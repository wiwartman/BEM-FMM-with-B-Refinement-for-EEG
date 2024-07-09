%   This is a mesh processor script: it computes necessary potential
%   integrals
%
%   Copyright SNM/WAW 2017-2024

% %%  Import geometry (exit if it exists)
% if exist('CombinedMeshP.mat', 'file')
%     tic
%     h = waitbar(0.5, 'Please wait - loading model data (if any)');
%     load CombinedMeshP.mat;
%     close(h);
%     disp([newline 'Model data loaded in ' num2str(toc) ' s']);
%     return;
% end

%%   Add accurate integration for electric field/electric potential on neighbor facets
%   Indexes into neighbor triangles
numThreads      = 4;       %   number of cores to be used
RnumberE        = 64;      %   number of neighbor triangles for analytical integration of electric field
RnumberP        = 64;      %   number of neighbor triangles for analytical integration of electric potential
ineighborE      = knnsearch(Center, Center, 'k', RnumberE);   % [1:N, 1:RnumberE]
ineighborP      = knnsearch(Center, Center, 'k', RnumberP);   % [1:N, 1:RnumberP]
ineighborE      = ineighborE';          %   do transpose  
ineighborP      = ineighborP';          %   do transpose  

tic
ppool           = gcp('nocreate');
if isempty(ppool)
    tic
    parpool(numThreads);
    disp([newline 'Started parallel pool in ' num2str(toc) ' s']);
end

%[EC, PC] = meshneighborints(P, t, normals, Area, Center, RnumberE, RnumberP, ineighborE, ineighborP, numThreads);
EC = meshneighborints_En(P, t, normals, Area, Center, RnumberE, ineighborE);
PC = meshneighborints_P(P, t, normals, Area, Center, RnumberP, ineighborP);


%%   Normalize sparse matrix EC by variable contrast (for speed up)
N   = size(Center, 1);
ii  = ineighborE;
jj  = repmat(1:N, RnumberE, 1); 
CO  = sparse(ii, jj, contrast(ineighborE));
EC  = CO.*EC;
