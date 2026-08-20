// Gmsh - Copyright (C) 1997-2026 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file in the Gmsh root directory for license information.
// Please report all issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#ifndef MESH_GREGION_MMG_H
#define MESH_GREGION_MMG_H

#include <vector>

class GRegion;

// Refine all of `regions` together in a single Mmg3d call, tagging each
// tetrahedron with its owning region so that interfaces between regions
// (which may share boundaries) are preserved as material interfaces by Mmg
// itself, instead of splitting the mesh into one independent per-region
// call beforehand.
void refineMeshMMG(std::vector<GRegion *> &regions);

#endif
