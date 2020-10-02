
function [A] = rotate_A_mex(A,angles,subMesh,A0)
% instead of recomputing all of A at each time / configuration, can simply rotate the entries
% corresponding to collocation points and elements on the same submesh, as
% long as the submesh appears to rigidly rotate over time

%leaves values alone for collocation point and element on different submeshes
%rotates values for collocation point and element on the input submesh to
%match current orientation angles

rotmat = rotation_matrix('z',angles(3)) * rotation_matrix('y',angles(2)) * rotation_matrix('x',angles(1));

blk = (kron(eye(subMesh.n_vert),rotmat));

A( (subMesh.global_indices.vert.start - 1)*3+1  :  subMesh.global_indices.vert.end * 3 , (subMesh.global_indices.vert.start - 1)*3+1  :  subMesh.global_indices.vert.end * 3) = blk * A0((subMesh.global_indices.vert.start - 1)*3+1  :  subMesh.global_indices.vert.end * 3 , (subMesh.global_indices.vert.start - 1)*3+1  :  subMesh.global_indices.vert.end * 3) * blk';
