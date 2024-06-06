%   This script plots the electric potential of the primary, secondary, or
%   the full field for any brain compartment surface
%
%   Copyright SNM/WAW 2017-2020

%%   Graphics
tissue_to_plot = 'Skin';
objectnumber    = find(strcmp(tissue, tissue_to_plot));
temp            = Ptot(Indicator==objectnumber);
PtotSkin        = temp;

%%  Digitize figure
figure;
step = 10;
temp = round(step*temp/max(temp)).*(max(temp))/step;
bemf2_graphics_surf_field(P, t, temp, Indicator, objectnumber);
title(strcat('Solution: Electric potential in V for:', tissue{objectnumber}));
view(-70, 70); colormap jet;