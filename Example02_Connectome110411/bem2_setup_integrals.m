%   This is a mesh processor script: it computes necessary potential
%   integrals
%
%   Copyright SNM/WAW 2017-2024

%%   Add accurate integration for electric field/electric potential on neighbor facets - global solution
%   Indexes into neighbor triangles
ineighborE      = knnsearch(Center, Center, 'k', RnumberE);   % [1:N, 1:RnumberE]
ineighborE      = ineighborE';          %   do transpose  
EC              = meshneighborints_En(P, t, normals, Area, Center, RnumberE, ineighborE);
tic
r1          = P(t(:, 1), :);
r2          = P(t(:, 2), :);
r3          = P(t(:, 3), :);
const       = 1/(4*pi); 
PC           = const*Area.*potint_self(r1, r2, r3); 
PotentialIntTime = toc


%%   Normalize sparse matrix EC by variable contrast (for speed up)
N   = size(Center, 1);
ii  = ineighborE;
jj  = repmat(1:N, RnumberE, 1); 
CO  = sparse(ii, jj, contrast(ineighborE));
EC  = CO.*EC;

%%   Local: add accurate integration for electric field on neighbor facets - local solution
%   Indexes into neighbor triangles
ineighborE      = knnsearch(L.Center, L.Center, 'k', L.RnumberE);   % [1:N, 1:RnumberE]
ineighborE      = ineighborE';          %   do transpose  

L.EC = meshneighborints_En(P, L.t, L.normals, L.Area, L.Center, L.RnumberE, ineighborE);

%%   Local: normalize sparse matrix L.EC by variable contrast (for speed up)
N   = size(L.Center, 1);
ii  = ineighborE;
jj  = repmat(1:N, L.RnumberE, 1); 
L.CO  = sparse(ii, jj, L.contrast(ineighborE));
L.EC  = L.CO.*L.EC;
