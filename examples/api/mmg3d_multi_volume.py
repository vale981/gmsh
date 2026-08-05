import gmsh
import sys

# Two unit cubes sharing a face, meshed with the MMG3D algorithm
# (Mesh.Algorithm3D = 7). This exercises MMG3D on a model with several
# volumes that share an interface, and checks that the resulting mesh is
# still conformal (no duplicated nodes on the shared face, no volume lost
# or gained).
#
# The mesh size is driven by an anisotropic metric (MathEvalAniso field),
# coarse in x/y and fine along z, so elements get stretched in x/y and
# refined in z. BAMG (Mesh.Algorithm = 7) is used for the 2D surface mesh
# since it is the only 2D algorithm that honors anisotropic metrics; MMG3D
# picks up the same background field for the 3D volume mesh.

gmsh.initialize()
gmsh.model.add("mmg3d_multi_volume")

gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1, 1)
gmsh.model.occ.addBox(1, 0, 0, 1, 1, 1, 2)
gmsh.model.occ.fragment([(3, 1)], [(3, 2)])
gmsh.model.occ.synchronize()

lc_xy = 0.3
lc_z = 0.1
gmsh.model.mesh.field.add("MathEvalAniso", 1)
gmsh.model.mesh.field.setString(1, "M11", f"1/{lc_xy}^2")
gmsh.model.mesh.field.setString(1, "M22", f"1/{lc_xy}^2")
gmsh.model.mesh.field.setString(1, "M33", f"1/{lc_z}^2")
gmsh.model.mesh.field.setString(1, "M12", "0")
gmsh.model.mesh.field.setString(1, "M13", "0")
gmsh.model.mesh.field.setString(1, "M23", "0")
gmsh.model.mesh.field.setAsBackgroundMesh(1)

# Let the field fully drive the size, as recommended when using a field
# with a large size gradient (see tutorial t10.py):
gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

gmsh.option.setNumber("Mesh.Algorithm", 7)  # BAMG, for the anisotropic 2D mesh
gmsh.option.setNumber("Mesh.Algorithm3D", 7)  # MMG3D

gmsh.model.mesh.generate(3)

# Sanity checks: the two volumes should still share the interface mesh
# (no duplicated coincident nodes), and the total tetrahedron volume
# should match the two unit cubes exactly.
node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
coords = node_coords.reshape((-1, 3))
seen = set()
duplicates = 0
for x, y, z in coords:
    key = (round(x, 9), round(y, 9), round(z, 9))
    if key in seen:
        duplicates += 1
    else:
        seen.add(key)
print("duplicated coincident nodes:", duplicates)

elem_types, _, elem_node_tags = gmsh.model.mesh.getElements(3)
vmap = {tag: coords[i] for i, tag in enumerate(node_tags)}
total_volume = 0.0
for etype, nodes in zip(elem_types, elem_node_tags):
    if etype != 4:  # 4-node tetrahedron
        continue
    nodes = nodes.reshape((-1, 4))
    for n in nodes:
        a, b, c, d = (vmap[t] for t in n)
        ab, ac, ad = b - a, c - a, d - a
        vol = abs(ab[0] * (ac[1] * ad[2] - ac[2] * ad[1]) -
                  ab[1] * (ac[0] * ad[2] - ac[2] * ad[0]) +
                  ab[2] * (ac[0] * ad[1] - ac[1] * ad[0])) / 6.0
        total_volume += vol
print("total tetrahedron volume (expected 2.0):", total_volume)

# Check that edges are indeed shorter along z than in the x/y plane.
edges = set()
for etype, nodes in zip(elem_types, elem_node_tags):
    if etype != 4:
        continue
    nodes = nodes.reshape((-1, 4))
    for n in nodes:
        for i in range(4):
            for j in range(i + 1, 4):
                edges.add((min(n[i], n[j]), max(n[i], n[j])))
dz, dxy = [], []
for a, b in edges:
    d = vmap[a] - vmap[b]
    dz.append(abs(d[2]))
    dxy.append((d[0]**2 + d[1]**2)**0.5)
print(f"average |dz| = {sum(dz) / len(dz):.4f}  "
      f"average |dxy| = {sum(dxy) / len(dxy):.4f} "
      f"(target ratio ~ {lc_z / lc_xy:.2f})")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
