[x, xi_eta, mesh_index, element_index, is_inside] = dist2mesh(Network.nodes, Mesh, [Metadata.geom.sphererad*1.1   max(sqrt(Mesh(1).area*4/sqrt(3)))*2; Metadata.geom.sphererad*1.1  max(sqrt(Mesh(2).area*4/sqrt(3)))*2  ]);


