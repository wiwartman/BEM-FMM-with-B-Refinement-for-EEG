%   This script plots the electric potential of the primary, secondary, or
%   the full field for any brain compartment surface
%
%   Copyright SNM/WAW 2017-2020

%%   Graphics
tissue_to_plot = 'Skin';
objectnumber    = find(strcmp(tissue, tissue_to_plot));
temp            = Ptot(Indicator==objectnumber);

figure;
bemf2_graphics_surf_field2(P, t, 1e6*temp, Indicator, objectnumber);
%title(strcat('Solution: Electric potential in uV for:', tissue{objectnumber}));
axis(1e-3*[-95 95 -95 95 -95 95]);
axis off