%   This script computes and plots the magnetic field magnitude
%   (or any of the components) for any brain compartment surface/interface
%    
%
%   Copyright SNM/WAW 2017-2024

%%   Compute the B-field for all surfaces/interfaces
tissue_to_plot = 'Skin';
objectnumber = find(strcmp(tissue, tissue_to_plot));    
Points = Center(Indicator==objectnumber, :); 
Normals = normals(Indicator==objectnumber, :); 

tic
d           = 10e-3; % Magnetometer distance from skin surface
obsPtsMag   = Points + d*Normals; % Observation points for magnetic field
[Idx, dist] = knnsearch(Center, obsPtsMag);
indexB      = dist<d/2;

difference      = condin - condout;
Bpri            = bemf3_inc_field_magnetic(strdipolemvector, strdipolemcenter, strdipolemstrength, obsPtsMag, mu0, prec);   
Bsec            = bemf5_volume_field_magnetic(obsPtsMag, Ptot, P, t, Center, Area, normals, difference, mu0, 0, 0, [], prec);
Bsec(indexB, :) = 0; %   Correction for very close points appearing due to deformations
Btotal          = Bpri + Bsec;

temp            = abs(sqrt(dot(Btotal, Btotal, 2)));
Btotal          = Bpri + Bsec;
temp            = abs(sqrt(dot(Btotal, Btotal, 2)));
BTime = toc

%%  Digitize figure
figure;
step = 15;
temp = round(step*temp/max(temp)).*(max(temp))/step;
bemf2_graphics_surf_field(P, t, temp, Indicator, objectnumber);
title(strcat('Solution: Magnetic field in T at 10 mm for:', tissue{objectnumber}));
view(-70, 70); colormap jet;