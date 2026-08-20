// Gmsh - Copyright (C) 1997-2026 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file in the Gmsh root directory for license information.
// Please report all issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#include <stdlib.h>
#include <vector>
#include "GmshConfig.h"
#include "GmshMessage.h"
#include "meshGRegion.h"
#include "meshGRegionHxt.h"
#include "meshGRegionNetgen.h"
#include "meshGRegionMMG.h"
#include "meshGFace.h"
#include "meshGFaceOptimize.h"
#include "meshGRegionBoundaryRecovery.h"
#include "meshGRegionDelaunayInsertion.h"
#include "meshRelocateVertex.h"
#include "meshUntangle.h"
#include "GModel.h"
#include "GRegion.h"
#include "GFace.h"
#include "GEdge.h"
#include "discreteFace.h"
#include "discreteEdge.h"
#include "MLine.h"
#include "MTriangle.h"
#include "MTetrahedron.h"
#include "MPyramid.h"
#include "MTrihedron.h"
#include "ExtrudeParams.h"
#include "OS.h"
#include "Context.h"
#include "SVector3.h"
#include <algorithm>
#include <cmath>

namespace {

// Average edge length of a tetrahedron, used as its local length scale.
double tetCharacteristicLength(MTetrahedron *t)
{
  double sum = 0;
  int n = 0;
  for(int i = 0; i < 4; i++) {
    for(int j = i + 1; j < 4; j++) {
      sum += t->getVertex(i)->distance(t->getVertex(j));
      n++;
    }
  }
  return n > 0 ? sum / n : 0;
}

// Repair near-zero-volume ("sliver") tetrahedra left over from boundary
// recovery by nudging one free (dim 3) vertex slightly off the plane of the
// other three. A handful of these can occur at numerically hard multi-way
// junctions between facets (more than two regions/facets meeting at a
// segment); Mmg3d flatly rejects any zero-volume tetrahedron as input
// rather than tolerating and later fixing it like it does other
// poor-quality elements, so such a tetrahedron must be repaired (or at
// least attempted) before reaching it.
//
// Tetrahedra with all 4 vertices constrained to a curve or surface (no
// free dim-3 vertex to move) are left untouched: repairing those would
// require moving a boundary vertex, which is out of scope here.
int repairSliverTetrahedra(GRegion *gr)
{
  // Reference length for detecting genuine degeneracy: the model's own
  // bounding box, not the tetrahedron's own edge lengths. An anisotropic
  // mesh legitimately contains many very flat tetrahedra (tiny volume
  // relative to their longest edge cubed), so comparing against the
  // tetrahedron's own size would flag those too; a truly degenerate
  // (near-coplanar) tetrahedron instead has a volume near floating-point
  // round-off relative to the coordinates involved, regardless of its
  // aspect ratio.
  SBoundingBox3d bbox = gr->bounds();
  double modelScale = bbox.diag();
  if(modelScale <= 0) modelScale = 1.;
  double volEps = 1.e-10 * modelScale * modelScale * modelScale;

  int nFixed = 0;
  for(MTetrahedron *t : gr->tetrahedra) {
    double L = tetCharacteristicLength(t);
    if(L <= 0) continue;
    double vol = std::fabs(t->getVolume());
    if(vol > volEps) continue; // not degenerate

    // Prefer the least-constrained vertex (highest dim() first): moving a
    // free interior (dim 3) point is exact, moving one on a curve/surface
    // is a tiny approximation (it ends up slightly off the exact CAD
    // geometry), but degenerate tetrahedra from segment/facet Steiner
    // point insertion typically have all 4 vertices boundary-constrained,
    // so requiring a dim-3 vertex would leave them all unrepaired.
    MVertex *victim = nullptr;
    int victimDim = -1;
    for(int i = 0; i < 4; i++) {
      MVertex *v = t->getVertex(i);
      int dim = v->onWhat() ? v->onWhat()->dim() : -1;
      if(dim > victimDim) {
        victim = v;
        victimDim = dim;
      }
    }
    if(!victim) continue;

    MVertex *a = nullptr, *b = nullptr, *c = nullptr;
    for(int i = 0; i < 4; i++) {
      MVertex *v = t->getVertex(i);
      if(v == victim) continue;
      if(!a) a = v;
      else if(!b) b = v;
      else c = v;
    }

    SVector3 ab(b->x() - a->x(), b->y() - a->y(), b->z() - a->z());
    SVector3 ac(c->x() - a->x(), c->y() - a->y(), c->z() - a->z());
    SVector3 normal = crossprod(ab, ac);
    if(normal.norm() <= 0) continue;
    normal.normalize();

    double eps = 1.e-4 * L;
    victim->setXYZ(victim->x() + eps * normal.x(),
                   victim->y() + eps * normal.y(),
                   victim->z() + eps * normal.z());
    nFixed++;
  }
  return nFixed;
}

} // namespace

void splitQuadRecovery::add(const MFace &f, MVertex *v, GFace *gf)
{
  _quad[f] = v;
  if(v) {
    MFace f0(f.getVertex(0), f.getVertex(1), v);
    MFace f1(f.getVertex(1), f.getVertex(2), v);
    MFace f2(f.getVertex(2), f.getVertex(3), v);
    MFace f3(f.getVertex(3), f.getVertex(0), v);
    _tri[f0] = gf;
    _tri[f1] = gf;
    _tri[f2] = gf;
    _tri[f3] = gf;
  }
  else {
    MTriangle t0(f.getVertex(0), f.getVertex(1), f.getVertex(2));
    MTriangle t1(f.getVertex(0), f.getVertex(2), f.getVertex(3));
    double qual01 = std::min(t0.gammaShapeMeasure(), t1.gammaShapeMeasure());
    MTriangle t2(f.getVertex(1), f.getVertex(2), f.getVertex(3));
    MTriangle t3(f.getVertex(0), f.getVertex(1), f.getVertex(3));
    double qual23 = std::min(t2.gammaShapeMeasure(), t3.gammaShapeMeasure());
    if (qual01 > qual23) {
      MFace f0(f.getVertex(0), f.getVertex(1), f.getVertex(2));
      MFace f1(f.getVertex(0), f.getVertex(2), f.getVertex(3));
      _tri[f0] = gf;
      _tri[f1] = gf;
    }
    else {
      MFace f0(f.getVertex(1), f.getVertex(2), f.getVertex(3));
      MFace f1(f.getVertex(0), f.getVertex(1), f.getVertex(3));
      _tri[f0] = gf;
      _tri[f1] = gf;
    }
  }
}

int splitQuadRecovery::buildPyramids(GModel *gm)
{
  if(_quad.empty()) return 0;

  Msg::Info("Generating pyramids for hybrid mesh...");
  int npyram = 0;
  int ntrihedra = 0;
  for(auto it = gm->firstRegion(); it != gm->lastRegion(); it++) {
    GRegion *gr = *it;
    if(gr->meshAttributes.method == MESH_TRANSFINITE) continue;
    if(gr->isFullyDiscrete()) {
      continue;
    }
    ExtrudeParams *ep = gr->meshAttributes.extrude;
    if(ep && ep->mesh.ExtrudeMesh && ep->geo.Mode == EXTRUDED_ENTITY) continue;

    std::vector<GFace *> faces = gr->faces();
    for(std::size_t i = 0; i < faces.size(); i++) {
      GFace *gf = faces[i];
      for(std::size_t j = 0; j < gf->quadrangles.size(); j++) {
        auto it2 = _quad.find(gf->quadrangles[j]->getFace(0));
        if(it2 != _quad.end()) {
          if(it2->second) {
            npyram++;
            gr->pyramids.push_back(new MPyramid(
              it2->first.getVertex(0), it2->first.getVertex(1),
              it2->first.getVertex(2), it2->first.getVertex(3), it2->second));
            gr->mesh_vertices.push_back(it2->second);
            if(it2->second->onWhat()->dim() == 3) {
              Msg::Error(
                "Pyramid top vertex already classified on volume %d (!= %d) - "
                "non-manifold quad boundaries not supported yet",
                it2->second->onWhat()->tag(), gr->tag());
            }
            else {
              it2->second->setEntity(gr);
            }
          }
          else {
            ntrihedra++;
            gr->trihedra.push_back(new MTrihedron(
              it2->first.getVertex(1), it2->first.getVertex(2),
              it2->first.getVertex(3), it2->first.getVertex(0)));
          }
        }
      }
    }
  }
  Msg::Info("Done generating %d pyramids and %d trihedra for hybrid mesh",
            npyram, ntrihedra);
  return npyram + ntrihedra;
}

static void _deleteUnusedVertices(GRegion *gr)
{
  std::set<MVertex *, MVertexPtrLessThan> allverts;
  for(std::size_t i = 0; i < gr->tetrahedra.size(); i++) {
    for(int j = 0; j < 4; j++) {
      if(gr->tetrahedra[i]->getVertex(j)->onWhat() == gr)
        allverts.insert(gr->tetrahedra[i]->getVertex(j));
    }
  }
  for(std::size_t i = 0; i < gr->mesh_vertices.size(); i++) {
    // FIXME: investigate crash on exit (e.g. t16.geo)
    // if(allverts.find(gr->mesh_vertices[i]) == allverts.end())
    //   delete gr->mesh_vertices[i];
  }
  gr->mesh_vertices.clear();
  gr->mesh_vertices.insert(gr->mesh_vertices.end(), allverts.begin(),
                           allverts.end());
}

void MeshDelaunayVolume(std::vector<GRegion *> &regions)
{
  if(regions.empty()) return;

  if(CTX::instance()->mesh.algo3d == ALGO_3D_HXT) {
    if(meshGRegionHxt(regions) != 0) { Msg::Error("HXT 3D mesh failed"); }
    return;
  }

  if(CTX::instance()->mesh.algo3d != ALGO_3D_RTREE &&
     CTX::instance()->mesh.algo3d != ALGO_3D_DELAUNAY &&
     CTX::instance()->mesh.algo3d != ALGO_3D_INITIAL_ONLY &&
     CTX::instance()->mesh.algo3d != ALGO_3D_MMG3D)
    return;

  GRegion *gr = regions[0];
  std::vector<GFace *> faces = gr->faces();

  std::set<GFace *, GEntityPtrLessThan> allFacesSet;
  for(std::size_t i = 0; i < regions.size(); i++) {
    std::vector<GFace *> const &f = regions[i]->faces();
    std::vector<GFace *> const &f_e = regions[i]->embeddedFaces();
    allFacesSet.insert(f.begin(), f.end());
    allFacesSet.insert(f_e.begin(), f_e.end());
  }

  // replace faces with compounds if elements from compound surface meshes are
  // not reclassified on the original surfaces
  if(CTX::instance()->mesh.compoundClassify == 0) {
    std::set<GFace *, GEntityPtrLessThan> comp;
    for(auto it = allFacesSet.begin(); it != allFacesSet.end(); it++) {
      GFace *gf = *it;
      if(!gf->compoundSurface)
        comp.insert(gf);
      else if(gf->compoundSurface)
        comp.insert(gf->compoundSurface);
    }
    allFacesSet = comp;
  }

  std::vector<GFace *> allFaces(allFacesSet.begin(), allFacesSet.end());
  gr->set(allFaces);

  std::set<GEdge *, GEntityPtrLessThan> allEmbEdgesSet;
  for(std::size_t i = 0; i < regions.size(); i++) {
    std::vector<GEdge *> const &e = regions[i]->embeddedEdges();
    allEmbEdgesSet.insert(e.begin(), e.end());
  }
  std::vector<GEdge *> allEmbEdges(allEmbEdgesSet.begin(),
                                   allEmbEdgesSet.end());
  std::vector<GEdge *> oldEmbEdges = gr->embeddedEdges();
  gr->embeddedEdges() = allEmbEdges;

  std::set<GVertex *> allEmbVerticesSet;
  for(std::size_t i = 0; i < regions.size(); i++) {
    std::vector<GVertex *> const &e = regions[i]->embeddedVertices();
    allEmbVerticesSet.insert(e.begin(), e.end());
  }
  std::vector<GVertex *> allEmbVertices(allEmbVerticesSet.begin(),
                                        allEmbVerticesSet.end());
  std::vector<GVertex *> oldEmbVertices = gr->embeddedVertices();
  gr->embeddedVertices() = allEmbVertices;

  splitQuadRecovery sqr(CTX::instance()->mesh.optimizePyramids >= -2);
  bool success = meshGRegionBoundaryRecovery(gr, &sqr);

  // sort triangles in all model faces in order to be able to search in vectors
  auto itf = allFaces.begin();
  while(itf != allFaces.end()) {
    std::sort((*itf)->triangles.begin(), (*itf)->triangles.end(),
              compareMTriangleLexicographic());
    ++itf;
  }

  // restore set of faces and embedded edges/vertices
  if(CTX::instance()->mesh.compoundClassify == 0) {
    std::set<GFace *, GEntityPtrLessThan> comp;
    for(std::size_t i = 0; i < faces.size(); i++) {
      GFace *gf = faces[i];
      if(!gf->compoundSurface)
        comp.insert(gf);
      else if(gf->compoundSurface)
        comp.insert(gf->compoundSurface);
    }
    std::vector<GFace *> lcomp(comp.begin(), comp.end());
    gr->set(lcomp);
  }
  else {
    gr->set(faces);
  }
  gr->embeddedEdges() = oldEmbEdges;
  gr->embeddedVertices() = oldEmbVertices;

  if(!success) return;

  // now do insertion of points
  if(CTX::instance()->mesh.algo3d == ALGO_3D_MMG3D) {
    // Classify the (possibly merged, if regions share a boundary) Delaunay
    // tetrahedra back onto their individual GRegion. maxIter=1 with a huge
    // radius target makes this call classify only, without inserting any
    // new point. All regions are then handed to Mmg3d together in a single
    // call (see refineMeshMMG), tagged with their region so that Mmg
    // preserves the interfaces between them as material boundaries; this
    // avoids relying on each region's boundary-recovered triangulation
    // being independently self-consistent, which can break down at
    // junctions shared by more than two regions.
    // Iterate: nudging a vertex shared by several tetrahedra to fix one
    // sliver can leave (or reveal) another degenerate one touching it, so
    // a single pass isn't always enough. Stop once a pass finds nothing
    // left to fix, or after a bounded number of rounds to avoid looping
    // forever on a case this simple repair can't actually resolve.
    int nSliversFixedTotal = 0;
    for(int round = 0; round < 10; round++) {
      int nSliversFixed = repairSliverTetrahedra(gr);
      if(nSliversFixed == 0) break;
      nSliversFixedTotal += nSliversFixed;
    }
    if(nSliversFixedTotal > 0) {
      Msg::Info("Repaired %d sliver tetrahedron%s left over from boundary "
                "recovery",
                nSliversFixedTotal, (nSliversFixedTotal > 1) ? "s" : "");
    }
    insertVerticesInRegion(gr, 1, 1.e300, true, &sqr);
    refineMeshMMG(regions);
  }
  else if(CTX::instance()->mesh.algo3d != ALGO_3D_INITIAL_ONLY &&
	  CTX::instance()->mesh.algo3d != ALGO_3D_RTREE) {
    insertVerticesInRegion(gr, CTX::instance()->mesh.maxIterDelaunay3D, 1.,
                           true, &sqr);
    for(auto gr : regions) _deleteUnusedVertices(gr);

    int nHybrid = sqr.buildPyramids(gr->model());
    if(nHybrid && sqr.doWeCreatePyramids()) {
      //      Msg::Info("Optimizing pyramids for hybrid mesh...");
      gr->model()->setAllVolumesPositive();
      RelocateVerticesOfPyramids(regions, 3);
      // RelocateVertices(regions, 3);
      //      Msg::Info("Done optimizing pyramids for hybrid mesh");
    }

    // test:
    // bool createBoundaryLayerOneLayer(GRegion *gr, std::vector<GFace *> &
    // bls); createBoundaryLayerOneLayer(gr, allFaces);
  }
}

void deMeshGRegion::operator()(GRegion *gr)
{
  if(gr->isFullyDiscrete()) return;
  gr->deleteMesh();
}

void meshGRegion::operator()(GRegion *gr)
{
  gr->model()->setCurrentMeshEntity(gr);

  if(gr->isFullyDiscrete()) return;
  if(gr->meshAttributes.method == MESH_NONE) return;
  if(CTX::instance()->mesh.meshOnlyVisible && !gr->getVisibility()) return;
  if(CTX::instance()->mesh.meshOnlyEmpty && gr->getNumMeshElements()) return;

  ExtrudeParams *ep = gr->meshAttributes.extrude;
  if(ep && ep->mesh.ExtrudeMesh) return;

  // destroy the mesh if it exists
  deMeshGRegion dem;
  dem(gr);

  if(MeshTransfiniteVolume(gr)) return;

  if(CTX::instance()->mesh.algo3d != ALGO_3D_FRONTAL) {
    delaunay.push_back(gr);
  }
  else if(CTX::instance()->mesh.algo3d == ALGO_3D_FRONTAL) {
    meshGRegionNetgen(gr);
  }
}

void untangleMeshGRegion::operator()(GRegion *gr, bool always)
{
  gr->model()->setCurrentMeshEntity(gr);

  if(!always && gr->isFullyDiscrete()) return;

  // don't optimize extruded meshes
  if(gr->meshAttributes.method == MESH_TRANSFINITE) return;
  ExtrudeParams *ep = gr->meshAttributes.extrude;
  if(ep && ep->mesh.ExtrudeMesh && ep->geo.Mode == EXTRUDED_ENTITY) return;

  Msg::Info("Untangling volume %d", gr->tag());
  untangleMesh(gr);
}


void optimizeMeshGRegion::operator()(GRegion *gr, bool always)
{
  gr->model()->setCurrentMeshEntity(gr);

  if(!always && gr->isFullyDiscrete()) return;

  // don't optimize extruded meshes
  if(gr->meshAttributes.method == MESH_TRANSFINITE) return;
  ExtrudeParams *ep = gr->meshAttributes.extrude;
  if(ep && ep->mesh.ExtrudeMesh && ep->geo.Mode == EXTRUDED_ENTITY) return;

  Msg::Info("Optimizing volume %d", gr->tag());
  optimizeMesh(gr, qmTetrahedron::QMTET_GAMMA);
}

bool buildFaceSearchStructure(GModel *model, fs_cont &search,
                              bool onlyTriangles)
{
  search.clear();

  std::set<GFace *> faces_to_consider;
  auto rit = model->firstRegion();
  while(rit != model->lastRegion()) {
    std::vector<GFace *> _faces = (*rit)->faces();
    faces_to_consider.insert(_faces.begin(), _faces.end());
    rit++;
  }

  auto fit = faces_to_consider.begin();
  while(fit != faces_to_consider.end()) {
    for(std::size_t i = 0; i < (*fit)->getNumMeshElements(); i++) {
      MFace ff = (*fit)->getMeshElement(i)->getFace(0);
      if(!onlyTriangles || ff.getNumVertices() == 3) search[ff] = *fit;
    }
    ++fit;
  }
  return true;
}

bool buildEdgeSearchStructure(GModel *model, es_cont &search)
{
  search.clear();

  auto eit = model->firstEdge();
  while(eit != model->lastEdge()) {
    for(std::size_t i = 0; i < (*eit)->lines.size(); i++) {
      MVertex *p1 = (*eit)->lines[i]->getVertex(0);
      MVertex *p2 = (*eit)->lines[i]->getVertex(1);
      MVertex *p = std::min(p1, p2);
      search.insert(std::pair<MVertex *, std::pair<MLine *, GEdge *> >(
        p, std::pair<MLine *, GEdge *>((*eit)->lines[i], *eit)));
    }
    ++eit;
  }
  return true;
}

GFace *findInFaceSearchStructure(MVertex *p1, MVertex *p2, MVertex *p3,
                                 const fs_cont &search)
{
  MFace ff(p1, p2, p3);
  auto it = search.find(ff);
  if(it == search.end()) return nullptr;
  return it->second;
}

GFace *findInFaceSearchStructure(const MFace &ff, const fs_cont &search)
{
  auto it = search.find(ff);
  if(it == search.end()) return nullptr;
  return it->second;
}

GEdge *findInEdgeSearchStructure(MVertex *p1, MVertex *p2,
                                 const es_cont &search)
{
  MVertex *p = std::min(p1, p2);

  for(auto it = search.lower_bound(p); it != search.upper_bound(p); ++it) {
    MLine *l = it->second.first;
    GEdge *ge = it->second.second;
    if((l->getVertex(0) == p1 || l->getVertex(0) == p2) &&
       (l->getVertex(1) == p1 || l->getVertex(1) == p2))
      return ge;
  }
  return nullptr;
}
