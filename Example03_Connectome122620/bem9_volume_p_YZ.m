%   This script accurately computes and displays electric potential sampled
%   on a cross-section (sagittal plane) via the FMM method with accurate
%   neighbor integration
%
%   Copyright SNM/WAW 2018-2024

%%  Load/prepare data
planeABCD = [1 0 0 -1e-3*Xm]; %Analytical equation of the observation point plane for neighbor search acceleration

%% Define observation points in the cross-section
Ms = 200;
y = linspace(ymin, ymax, Ms);
z = linspace(zmin, zmax, Ms);
[Y0, Z0]  = meshgrid(y, z);
clear pointsYZ;
pointsYZ(:, 1) = X*ones(1, Ms^2);
pointsYZ(:, 2) = reshape(Y0, 1, Ms^2);
pointsYZ(:, 3) = reshape(Z0, 1, Ms^2);

%% Find the potential at each observation point
tic
pointsYZ       = pointsYZ;     % Convert back to m
Ppri           = zeros(Ms*Ms, 1);
Psec           = zeros(Ms*Ms, 1);
flag = 1;       %   precise integration
[~, Ppri]      = bemf3_inc_field_electric_plain_dipoles(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, 1e-3*pointsYZ, 1e-2*prec);
R = 2;  %   precise integration
Psec           = bemf5_volume_field_potential(1e-3*pointsYZ, c, P, t, Center, Area, normals, R, planeABCD, prec);
Ptotal         = Ppri + Psec;   
disp([newline 'Potential calculated in ' num2str(toc) ' s']);

%% Plot the potential in the cross-section
figure;
% Potential contour plot
temp      = 1e6*Ptotal;
th1 =  +2e-3*max(+temp);             %   in V
th2 =  -2e-3*max(-temp);             %   in V
levels      = 30;
bemf2_graphics_vol_field_log(temp, th1, th2, levels, y, z);
xlabel('Distance y, mm');
ylabel('Distance z, mm');
title(strcat('Potential (uV), ', label, '-in the sagittal plane'));

% Tissue boundaries
for m = countYZ
    edges           = EofYZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofYZ{m}(:, 2);       %   this is for the contour  
    points(:, 2)    = +PofYZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.5);    %   this is contour plot
end

% Dipoles
hold on
bemf1_graphics_dipole(1e3*strdipolePplus, 1e3*strdipolePminus, strdipoleCurrent, 3);

% General settings 
axis 'equal';  axis 'tight';     
colormap jet;
axis([ymin ymax zmin zmax]);
grid on; set(gcf,'Color','White');