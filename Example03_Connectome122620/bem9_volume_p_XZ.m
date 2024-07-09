%   This script accurately computes and displays electric potential sampled
%   on a cross-section (coronal plane) via the FMM method with accurate
%   neighbor integration
%
%   Copyright SNM/WAW 2018-2024

%%  Load/prepare data
planeABCD = [0 1 0 -1e-3*Ym]; %Analytical equation of the observation point plane for neighbor search acceleration

%% Define observation points in the cross-section
Ms = 300;
x = linspace(xmin, xmax, Ms);
z = linspace(zmin, zmax, Ms);
[X0, Z0]  = meshgrid(x, z);
clear pointsXZ;
pointsXZ(:, 1) = reshape(X0, 1, Ms^2);
pointsXZ(:, 2) = Y*ones(1, Ms^2);
pointsXZ(:, 3) = reshape(Z0, 1, Ms^2);  

%% Find the potential at each observation point
tic
pointsXZ       = pointsXZ;     % Convert back to m
Ppri           = zeros(Ms*Ms, 1);
Psec           = zeros(Ms*Ms, 1);
[~, Ppri]      = bemf3_inc_field_electric_plain_dipoles(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, 1e-3*pointsXZ, 1e-2*prec);
R = 2;  %   precise integration
Psec           = bemf5_volume_field_potential(1e-3*pointsXZ, c, P, t, Center, Area, normals, R, planeABCD, prec);
Ptotal         = Ppri + Psec;   
disp([newline 'Potential calculated in ' num2str(toc) ' s']);

%% Plot potential in the cross-section
figure;
% Potential contour plot
temp      = 1e6*Ptotal;
th1 =  +2e-3*max(+temp);             %   in V
th2 =  -2e-3*max(-temp);             %   in V
levels      = 30;
bemf2_graphics_vol_field_log(temp, th1, th2, levels, x, z);
xlabel('Distance x, mm');
ylabel('Distance z, mm');
title(strcat('Potential (uV), ', label, '-in the coronal plane'));
 
% Tissue boundaries
for m = countXZ
    edges           = EofXZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofXZ{m}(:, 1);       %   this is for the contour  
    points(:, 2)    = +PofXZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.5);    %   this is contour plot
end

% Dipoles
hold on
bemf1_graphics_dipole(1e3*strdipolePplus, 1e3*strdipolePminus, strdipoleCurrent, 2);

% General settings 
axis 'equal';  axis 'tight';     
colormap jet;
axis([xmin xmax zmin zmax]);
grid on; set(gcf,'Color','White');