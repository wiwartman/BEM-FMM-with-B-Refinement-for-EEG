%   This script plots the induced surface charge density for
%   any brain compartment surface (plots the density)
%
%   Copyright SNM/WAW 2018-2020

%%   Graphics
tissue_to_plot = 'Skin';
objectnumber    = find(strcmp(tissue, tissue_to_plot));
temp            = eps0*c(Indicator==objectnumber);  % the real charge density is eps0*c

scale = 0.5*max(abs(temp));
temp(temp>+scale) = +scale;
temp(temp<-scale) = -scale;

figure;
bemf2_graphics_surf_field(P, t, temp, Indicator, objectnumber);
%title(strcat('Solution: Surface charge density in C/m^2 for:', tissue{objectnumber}));
axis off;