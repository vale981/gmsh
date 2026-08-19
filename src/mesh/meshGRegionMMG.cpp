// Gmsh - Copyright (C) 1997-2026 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file in the Gmsh root directory for license information.
// Please report all issues on https://gitlab.onelab.info/gmsh/gmsh/issues.
//
// Contributor(s):
//   Algiane Froehly
//

#include "GmshConfig.h"
#include "GmshMessage.h"
#include "meshGRegionMMG.h"

#if defined(HAVE_MMG)

#include <algorithm>
#include <map>
#include <set>
#include <vector>
#include "GEntity.h"
#include "GRegion.h"
#include "GFace.h"
#include "MTetrahedron.h"
#include "MTriangle.h"
#include "MVertex.h"
#include "BackgroundMeshTools.h"
#include "Context.h"

extern "C" {
#include <mmg/libmmg.h>
}

namespace {

  // The union of the boundary faces of all `regions`, without double-counting
  // a face shared between two regions (an internal interface).
  std::vector<GFace *> boundaryFaces(const std::vector<GRegion *> &regions)
  {
    std::set<GFace *, GEntityPtrLessThan> facesSet;
    for(GRegion *gr : regions) {
      std::vector<GFace *> const &f = gr->faces();
      facesSet.insert(f.begin(), f.end());
    }
    return std::vector<GFace *>(facesSet.begin(), facesSet.end());
  }

} // namespace

static void MMG2gmsh(std::vector<GRegion *> &regions, MMG5_pMesh mmg,
                     std::map<int, MVertex *> &mmg2gmsh)
{
  std::map<int, GRegion *> tagToRegion;
  for(GRegion *gr : regions) tagToRegion[gr->tag()] = gr;

  int np, ne, nt, na, ref;
  double cx, cy, cz;

  if(MMG3D_Get_meshSize(mmg, &np, &ne, nullptr, &nt, nullptr, &na) != 1)
    Msg::Error("Mmg3d: unable to get mesh size");

  // Store the nodes from the Mmg structures into gmsh vertices. A vertex
  // that already existed before remeshing (found in mmg2gmsh) is reused as
  // is; a brand-new (Steiner) vertex is created without an owning entity
  // yet, since -- now that a single combined mesh can span several
  // GRegions -- its owning region is only known once a tetrahedron using it
  // is seen below.
  //
  // TODO: when MMG is allowed to modify the surface mesh, reclassify nodes
  // accordingly
  std::vector<MVertex *> kToMVertex(np + 1, nullptr);
  std::set<MVertex *> unattached;
  for(int k = 1; k <= np; k++) {
    if(MMG3D_Get_vertex(mmg, &cx, &cy, &cz, &ref, nullptr, nullptr) != 1)
      Msg::Error("Mmg3d: unable to get vertex %d", k);

    auto it = mmg2gmsh.find(ref);
    if(it == mmg2gmsh.end()) {
      MVertex *v = new MVertex(cx, cy, cz, nullptr);
      kToMVertex[k] = v;
      unattached.insert(v);
    }
    else {
      kToMVertex[k] = it->second;
    }
  }

  // Store the tets from the Mmg structures into the Gmsh structures,
  // routing each one to the GRegion matching its reference, and attaching
  // any brand-new vertex to the first region found using it.
  for(int k = 1; k <= ne; k++) {
    int v1mmg, v2mmg, v3mmg, v4mmg;
    if(MMG3D_Get_tetrahedron(mmg, &v1mmg, &v2mmg, &v3mmg, &v4mmg, &ref,
                             nullptr) != 1)
      Msg::Error("Mmg3d: unable to get tetrahedron %d", k);

    MVertex *v1 = kToMVertex[v1mmg];
    MVertex *v2 = kToMVertex[v2mmg];
    MVertex *v3 = kToMVertex[v3mmg];
    MVertex *v4 = kToMVertex[v4mmg];
    if(!v1 || !v2 || !v3 || !v4) {
      Msg::Error("Mmg3d: unknown vertex in tetrahedron %d", k);
      continue;
    }

    auto rit = tagToRegion.find(ref);
    if(rit == tagToRegion.end()) {
      Msg::Error("Mmg3d: tetrahedron %d has unknown region reference %d", k,
                 ref);
      continue;
    }
    GRegion *gr = rit->second;
    gr->tetrahedra.push_back(new MTetrahedron(v1, v2, v3, v4));

    for(MVertex *v : {v1, v2, v3, v4}) {
      if(unattached.erase(v)) {
        gr->mesh_vertices.push_back(v);
        v->setEntity(gr);
      }
    }
  }

#if 0
  // Store the triangles from the Mmg structures into the Gmsh structures

  // TODO: allow MMG to remesh (some of the) input surfaces; then store the new
  // mesh
  for(int k = 1; k <= nt; k++) {
    int v1mmg, v2mmg, v3mmg;
    if(MMG3D_Get_triangle(mmg, &v1mmg, &v2mmg, &v3mmg, &ref, nullptr) != 1)
      Msg::Error("Mmg3d: unable to get triangle %d", k);

    MVertex *v1 = kToMVertex[v1mmg];
    MVertex *v2 = kToMVertex[v2mmg];
    MVertex *v3 = kToMVertex[v3mmg];
    if(!v1 || !v2 || !v3) {
      Msg::Error("Mmg3d: unknown vertex in triangle %d", k);
    }
    else{
      // ...
    }
  }
#endif
}

// Returns false if any Mmg "Set" call failed, meaning `mmg`/`sol` were left
// partially initialized and must not be handed to MMG3D_mmg3dlib (doing so
// can crash deep inside Mmg instead of failing cleanly).
static bool gmsh2MMG(std::vector<GRegion *> &regions, MMG5_pMesh mmg,
                     MMG5_pSol sol, std::map<int, MVertex *> &mmg2gmsh)
{
  bool ok = true;

  // Count mesh vertices across all regions
  std::set<MVertex *> allVertices;
  for(GRegion *gr : regions) {
    for(std::size_t i = 0; i < gr->tetrahedra.size(); i++) {
      allVertices.insert(gr->tetrahedra[i]->getVertex(0));
      allVertices.insert(gr->tetrahedra[i]->getVertex(1));
      allVertices.insert(gr->tetrahedra[i]->getVertex(2));
      allVertices.insert(gr->tetrahedra[i]->getVertex(3));
    }
  }
  int np = allVertices.size();

  // Boundary triangles: the union of all regions' bounding faces, so a
  // face shared between two regions (an internal interface) is only set
  // once; the differing tetrahedron references on either side (set below)
  // are enough for Mmg to recognize and preserve it as a material
  // interface on its own.
  std::vector<GFace *> f = boundaryFaces(regions);
  int nt = 0;
  for(auto it = f.begin(); it != f.end(); ++it) nt += (*it)->triangles.size();

  // TODO: also import mesh lines

  // Get mesh tetrahedra
  int ne = 0;
  for(GRegion *gr : regions) ne += gr->tetrahedra.size();

  if(MMG3D_Set_meshSize(mmg, np, ne, 0, nt, 0, 0) != 1) {
    Msg::Error("Mmg3d: unable to set mesh size");
    return false;
  }

  if(MMG3D_Set_solSize(mmg, sol, MMG5_Vertex, np, MMG5_Tensor) != 1) {
    Msg::Error("Mmg3d: unable to set metric size");
    return false;
  }

  std::map<MVertex *, std::pair<double, int>> LCS;
  for(auto it = f.begin(); it != f.end(); ++it) {
    for(unsigned int i = 0; i < (*it)->triangles.size(); i++) {
      MTriangle *t = (*it)->triangles[i];
      double L = t->maxEdge();
      for(int k = 0; k < 3; k++) {
        MVertex *v = t->getVertex(k);
        auto itv = LCS.find(v);
        if(itv != LCS.end()) {
          itv->second.first += L;
          itv->second.second++;
        }
        else {
          LCS[v] = std::make_pair(L, 1);
        }
      }
    }
  }

  int k = 1;
  std::map<int, int> gmsh2mmg_num;
  for(auto it = allVertices.begin(); it != allVertices.end(); ++it) {
    if(MMG3D_Set_vertex(mmg, (*it)->x(), (*it)->y(), (*it)->z(),
                        (*it)->getNum(), k) != 1) {
      Msg::Error("Mmg3d: unable to set vertex %d", k);
      ok = false;
    }

    gmsh2mmg_num[(*it)->getNum()] = k;

    MVertex *v = *it;
    double U = 0, V = 0;
    if(!v->onWhat()) continue;

    if(v->onWhat()->dim() == 1) { v->getParameter(0, U); }
    else if(v->onWhat()->dim() == 2) {
      v->getParameter(0, U);
      v->getParameter(1, V);
    }

    // double lc = BGM_MeshSize(v->onWhat(), U, V, v->x(), v->y(), v->z());
    SMetric3 m = BGM_MeshMetric(v->onWhat(), U, V, v->x(), v->y(), v->z());

    auto itv = LCS.find(v);
    if(itv != LCS.end()) {
      mmg2gmsh[(*it)->getNum()] = *it;
      // if (Extend2dMeshIn3dVolumes()){
      double LL = itv->second.first / itv->second.second;
      SMetric3 l4(1. / (LL * LL));
      SMetric3 MM = intersection_conserve_mostaniso(l4, m);
      m = MM;
      // lc = std::min(LL,lc);
      // }
    }

    if(MMG3D_Set_tensorSol(sol, m(0, 0), m(1, 0), m(2, 0), m(1, 1), m(2, 1),
                           m(2, 2), k) != 1) {
      Msg::Error("Mmg3d: unable to set solution %d", k);
      ok = false;
    }
    k++;
  }

  k = 1;
  for(GRegion *gr : regions) {
    for(std::size_t i = 0; i < gr->tetrahedra.size(); i++) {
      if(MMG3D_Set_tetrahedron(
           mmg, gmsh2mmg_num[gr->tetrahedra[i]->getVertex(0)->getNum()],
           gmsh2mmg_num[gr->tetrahedra[i]->getVertex(1)->getNum()],
           gmsh2mmg_num[gr->tetrahedra[i]->getVertex(2)->getNum()],
           gmsh2mmg_num[gr->tetrahedra[i]->getVertex(3)->getNum()], gr->tag(),
           k) != 1) {
        Msg::Error("Mmg3d: unable to set tetrahedron %d (region %d): "
                   "degenerate (zero-volume) tetrahedron in the classified "
                   "mesh",
                   k, gr->tag());
        ok = false;
      }
      k++;
    }
  }

  k = 1;
  for(auto it = f.begin(); it != f.end(); ++it) {
    for(unsigned int i = 0; i < (*it)->triangles.size(); i++) {
      if(MMG3D_Set_triangle(
           mmg, gmsh2mmg_num[(*it)->triangles[i]->getVertex(0)->getNum()],
           gmsh2mmg_num[(*it)->triangles[i]->getVertex(1)->getNum()],
           gmsh2mmg_num[(*it)->triangles[i]->getVertex(2)->getNum()],
           (*it)->tag(), k) != 1) {
        Msg::Error("Mmg3d: unable to set triangle %d", k);
        ok = false;
      }
      k++;
    }
  }

  return ok;
}

static void updateSizes(std::vector<GRegion *> &regions, MMG5_pMesh mmg,
                        MMG5_pSol sol, std::map<int, MVertex *> &mmg2gmsh)
{
  std::map<int, GRegion *> tagToRegion;
  for(GRegion *gr : regions) tagToRegion[gr->tag()] = gr;

  std::vector<GFace *> f = boundaryFaces(regions);

  std::map<MVertex *, std::pair<double, int>> LCS;
  // if (Extend2dMeshIn3dVolumes()){
  for(auto it = f.begin(); it != f.end(); ++it) {
    for(unsigned int i = 0; i < (*it)->triangles.size(); i++) {
      MTriangle *t = (*it)->triangles[i];
      double L = t->maxEdge();
      for(int k = 0; k < 3; k++) {
        MVertex *v = t->getVertex(k);
        auto itv = LCS.find(v);
        if(itv != LCS.end()) {
          itv->second.first += L;
          itv->second.second++;
        }
        else {
          LCS[v] = std::make_pair(L, 1);
        }
      }
    }
  }
  // }

  int np, ne;
  MMG3D_Get_meshSize(mmg, &np, &ne, nullptr, nullptr, nullptr, nullptr);

  // Determine, for every brand-new (interior) vertex with no pre-existing
  // entity, which GRegion owns it via any tetrahedron using it: Mmg
  // preserves the boundary between differently-referenced tetrahedra (see
  // gmsh2MMG), so every tetrahedron touching a genuinely interior vertex
  // shares the same reference, and any one of them suffices to identify the
  // owning region.
  std::vector<GRegion *> vertexRegion(np + 1, nullptr);
  for(int k = 1; k <= ne; k++) {
    int v1, v2, v3, v4, ref;
    if(MMG3D_Get_tetrahedron(mmg, &v1, &v2, &v3, &v4, &ref, nullptr) != 1) {
      // v1..v4 are left uninitialized on failure: skip rather than use them
      // as (garbage) indices into vertexRegion below.
      Msg::Error("Mmg3d: unable to get tetrahedron %d", k);
      continue;
    }
    auto rit = tagToRegion.find(ref);
    GRegion *owner = (rit != tagToRegion.end()) ? rit->second : nullptr;
    for(int vv : {v1, v2, v3, v4}) {
      if(vv >= 1 && vv <= np && !vertexRegion[vv]) vertexRegion[vv] = owner;
    }
  }

  for(int k = 1; k <= np; k++) {
    double cx, cy, cz;
    if(MMG3D_Get_vertex(mmg, &cx, &cy, &cz, nullptr, nullptr, nullptr) != 1)
      Msg::Error("Mmg3d: unable to get vertex %d", k);

    auto it = mmg2gmsh.find(k);
    GEntity *ge = (it != mmg2gmsh.end() && it->second->onWhat()) ?
                    it->second->onWhat() :
                    static_cast<GEntity *>(vertexRegion[k]);
    if(!ge) continue;

    SMetric3 m = BGM_MeshMetric(ge, 0, 0, cx, cy, cz);

    if(it != mmg2gmsh.end() && Extend2dMeshIn3dVolumes()) {
      auto itv = LCS.find(it->second);
      if(itv != LCS.end()) {
        double LL = itv->second.first / itv->second.second;
        // printf("adding a size %g\n",LL);
        SMetric3 l4(1. / (LL * LL));
        SMetric3 MM = intersection_conserve_mostaniso(l4, m);
        m = MM;
      }
    }
    if(m.determinant() < 1.e-30) {
      m(0, 0) += 1.e-12;
      m(1, 1) += 1.e-12;
      m(2, 2) += 1.e-12;
    }

    if(MMG3D_Set_tensorSol(sol, m(0, 0), m(1, 0), m(2, 0), m(1, 1), m(2, 1),
                           m(2, 2), k) != 1)
      Msg::Error("Mmg3d: unable to set solution %d", k);
  }
}

void refineMeshMMG(std::vector<GRegion *> &regions)
{
  if(regions.empty()) return;

  MMG5_pMesh mmg = nullptr;
  MMG5_pSol sol = nullptr;

  std::map<int, MVertex *> mmg2gmsh;

  // Mmg structures allocations
  MMG3D_Init_mesh(MMG5_ARG_start, MMG5_ARG_ppMesh, &mmg, MMG5_ARG_ppMet, &sol,
                  MMG5_ARG_end);

  // Store the Gmsh mesh (all regions at once) into the Mmg structures
  if(!gmsh2MMG(regions, mmg, sol, mmg2gmsh)) {
    Msg::Error("Mmg3d: failed to build the input mesh (see above); leaving "
               "the classified mesh unrefined instead of handing a "
               "partially-built mesh to Mmg3d");
    MMG3D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mmg, MMG5_ARG_ppMet, &sol,
                   MMG5_ARG_end);
    return;
  }

  int iterMax = 10;

// #define DEBUG
#ifdef DEBUG
  char test0[] = "init.mesh";
  MMG3D_saveMesh(mmg, test0);
  MMG3D_saveSol(mmg, sol, test0);
#endif

  for(int ITER = 0; ITER < iterMax; ITER++) {
    int nT, nTnow, np;

    MMG3D_Get_meshSize(mmg, nullptr, &nT, nullptr, nullptr, nullptr, nullptr);

    // Mmg parameters : verbosity + nosurf option
    int verb_mmg = (Msg::GetVerbosity() < 2) ? -1 : Msg::GetVerbosity() - 4;
    if(MMG3D_Set_iparameter(mmg, sol, MMG3D_IPARAM_verbose, verb_mmg) != 1)
      Msg::Error("Mmsg3d: unable to set verbosity");

    // Set the nosurf parameter to 1 to preserve the boundaries. This also
    // preserves the interfaces between differently-referenced regions: by
    // default (opnbdy off), Mmg already treats any triangle between two
    // tetrahedra of different references as a boundary to keep, so a
    // multi-region mesh remeshed in one combined call still has its
    // internal material interfaces respected.
    if(MMG3D_Set_iparameter(mmg, sol, MMG3D_IPARAM_nosurf, 1) != 1)
      Msg::Error("Mmg3d: unable to preserve the boundaries");

    // Set the hausdorff parameter, scaled by the largest region's bounding
    // box (matches the previous per-region behaviour when there is only
    // one region).
    double sqrt3Inv = 0.57735026919;
    double diag = 0;
    for(GRegion *gr : regions) diag = std::max(diag, gr->bounds().diag());
    double hausd = 0.01 * sqrt3Inv * diag;

    if(MMG3D_Set_dparameter(mmg, sol, MMG3D_DPARAM_hausd, hausd) != 1) {
      Msg::Error("Mmg3d: unable to set the hausdorff parameter");
    }

    if(MMG3D_mmg3dlib(mmg, sol) != MMG5_SUCCESS) {
      Msg::Error("Mmg3d: failed (iteration %d)", ITER);
    }
    else {
      MMG3D_Get_meshSize(mmg, &np, &nTnow, nullptr, nullptr, nullptr, nullptr);
      Msg::Info("Mmg3d: success (iteration %d) - %d nodes %d tetrahedra", ITER,
                np, nTnow);

      // Here we should interact with BGM
      updateSizes(regions, mmg, sol, mmg2gmsh);

      if(fabs((double)(nTnow - nT)) < 0.05 * nT) break;
    }
  }

#ifdef DEBUG
  char test[] = "end.mesh";
  MMG3D_saveMesh(mmg, test);
  MMG3D_saveSol(mmg, sol, test);
#endif

  for(GRegion *gr : regions) {
    gr->deleteVertexArrays();
    for(unsigned int i = 0; i < gr->tetrahedra.size(); ++i)
      delete gr->tetrahedra[i];
    gr->tetrahedra.clear();
    for(unsigned int i = 0; i < gr->mesh_vertices.size(); ++i)
      delete gr->mesh_vertices[i];
    gr->mesh_vertices.clear();
  }

  // Store the Mmg mesh into the Gmsh structures, routing each tetrahedron
  // and vertex back to its owning region
  MMG2gmsh(regions, mmg, mmg2gmsh);

  // Free the Mmg structure
  MMG3D_Free_all(MMG5_ARG_start, MMG5_ARG_ppMesh, &mmg, MMG5_ARG_ppMet, &sol,
                 MMG5_ARG_end);
}

#else

void refineMeshMMG(std::vector<GRegion *> &regions)
{
  Msg::Warning("This version of Gmsh is not compiled with MMG support: "
               "skipping refinement");
}

#endif
