
function P = elastic_getP_new(n,npmlz,npmlx,Iz,Ix)
% sampling operator, selects gridpoints from 2D grid:
%
% use:
%   P = getP(n,Iz,Ix);
%
% input:
%   n     - [nz,nx] number of gridpoints in z and x direction
%   Iz,Ix - indices
%
% ouput:
%   P - sparse matrix
%
Nx = n(2) + 2*npmlx;
Nz = n(1) + 2*npmlz;
I1 = speye(Nz);
I2 = speye(Nx);
P  = kron(I1(npmlz+Iz,:),I2(npmlx+Ix,:));

