
function [A] = rotate_A0(A,angles,subMesh,A0)
% instead of recomputing all of A at each time / configuration, can simply rotate the entries
% corresponding to collocation points and elements on the same submesh, as
% long as the submesh appears to rigidly rotate over time

%leaves values alone for collocation point and element on different submeshes
%rotates values for collocation point and element on the input submesh to
%match current orientation angles

%algorithm to rotate a tensor stolen from some internet forum.  will need to find again at
%some point....

rotmat = rotation_matrix('z',angles(3)) * rotation_matrix('y',angles(2)) * rotation_matrix('x',angles(1));

blk = kron(speye(subMesh.n_vert),rotmat);  %using sparse is *a lot* faster than not - I checked.  also, mexed version of this function is slower - turns out shatlab is actually pretty good at some things

A( (subMesh.indices.glob.unq_bounds.vert(1) - 1)*3+1  :  subMesh.indices.glob.unq_bounds.vert(2) * 3 , ...
   (subMesh.indices.glob.unq_bounds.vert(1) - 1)*3+1  :  subMesh.indices.glob.unq_bounds.vert(2) * 3) ...
   = blk *            A0((subMesh.indices.glob.unq_bounds.vert(1) - 1)*3+1  :  subMesh.indices.glob.unq_bounds.vert(2) * 3 , ...
                         (subMesh.indices.glob.unq_bounds.vert(1) - 1)*3+1  :  subMesh.indices.glob.unq_bounds.vert(2) * 3)        * blk';



