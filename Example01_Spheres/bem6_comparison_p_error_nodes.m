%   This script compares analytical and numerical solutions for the
%   potential.
%   An infinitesimally short dipole within a four-layer sphere is considered
%
%   Copyright SNM 2018-2020

radfactor     = [0.92 0.86 0.80 0.78];          %     skin skull csf gm 

%   Facets
PotNum = Ptot;         
tissue_to_plot  = 'Skin';
objectnumber    = find(strcmp(tissue, tissue_to_plot));
index           = Indicator==objectnumber;
PotNum          = PotNum(index);

%   Nodes
PSKIN   = P(1:max(max(t(index, :))), :);
tSKIN   = t(index, :);
DT      = triangulation(tSKIN, PSKIN); 
V       = vertexAttachments(DT);
PotVNum = zeros(size(PSKIN, 1), 1);
for m = 1:size(PSKIN, 1)
    PotVNum(m) = mean(PotNum(V{m}));     
end

%   Analytical solution
tic
M       = size(strdipolePplus, 1);
PotAnl  = zeros(size(PSKIN, 1), 1);

%parpool(20);
for m = 1:M
    temp    = a_p_4layer_infinite(strdipoleCurrent(m), strdipolePplus(m, :), strdipolePminus(m, :), radfactor, cond(1:4), PSKIN);
    temp    = real(temp);
    PotAnl  = PotAnl + temp;
end
%delete(gcp('nocreate'));
toc

f1 = figure;
bemf2_graphics_surf_field(P, t, PotNum, Indicator, objectnumber);
title('numerical')

f2 = figure;
title('analytical')
PotAnlt     = 1/3*(PotAnl(tSKIN(:, 1), :) + PotAnl(tSKIN(:, 2), :) + PotAnl(tSKIN(:, 3), :));
bemf2_graphics_surf_field(P, t, PotAnlt, Indicator, objectnumber);

Error_2norm = norm(PotVNum - PotAnl)/norm(PotAnl)
Error_RDM   = 0.5*norm(PotVNum/norm(PotVNum) - PotAnl/norm(PotAnl))
