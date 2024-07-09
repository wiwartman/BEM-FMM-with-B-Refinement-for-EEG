%   This script accurately computes and displays electric fields sampled on
%   a cross-section (transverse plane) via the FMM method with accurate
%   neighbor integration
%
%   Copyright SNM/WAW 2017-2024

%%  Load/prepare data
planeABCD = [0 0 1 -1e-3*Zm]; %Analytical equation of the observation point plane for neighbor search acceleration

%%  Post processing parameters
component   = 4;        %   field component to be plotted (1, 2, 3 or x, y, z, or 4 - total) 
temp        = ['x' 'y' 'z' 't'];
label       = temp(component);

%%  Define observation points in the cross-section
Ms = 200;
x = linspace(xmin, xmax, Ms);
y = linspace(ymin, ymax, Ms);
[X0, Y0]  = meshgrid(x, y);
clear pointsXY;
pointsXY(:, 1) = reshape(X0, 1, Ms^2);
pointsXY(:, 2) = reshape(Y0, 1, Ms^2);  
pointsXY(:, 3) = Zm*ones(1, Ms^2);

%%  Find the E-field at each observation point in the cross-section
tic
pointsXY       = pointsXY;     % Convert back to m
Epri           = zeros(Ms*Ms, 3);
Esec           = zeros(Ms*Ms, 3);
[Epri, ~]      = bemf3_inc_field_electric_plain_dipoles(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, 1e-3*pointsXY, 1e-2*prec);
R = 2;  %   precise integration
Esec           = bemf5_volume_field_electric(1e-3*pointsXY, c, P, t, Center, Area, normals, R, planeABCD, prec);
Etotal         = Epri + Esec; 
disp([newline 'E-field calculated in ' num2str(toc) ' s']);


%%  Plot the E-field in the cross-section
figure;
%  E-field contour plot
if component == 4
    temp      = abs(sqrt(dot(Etotal, Etotal, 2)));
else
    temp      = abs(Etotal(:, component));
end
th1 = +1e-3*max(temp);             %   in V/m
th2 =  min(temp);                  %   in V/m
levels      = 10;
bemf2_graphics_vol_field_log(temp, th1, th2, levels, x, y);
xlabel('Distance x, mm');
ylabel('Distance y, mm');
title(strcat('E-field (V/m), ', label, '-component in the transverse plane.'));

% Tissue boundaries
for m = countXY
    edges           = EofXY{m};             %   this is for the contour
    points          = [];
    points(:, 1)    = +PofXY{m}(:, 1);       %   this is for the contour  
    points(:, 2)    = +PofXY{m}(:, 2);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.5);    %   this is contour plot
end

% Dipoles
hold on
bemf1_graphics_dipole(1e3*strdipolePplus, 1e3*strdipolePminus, strdipoleCurrent, 1);

% General settings 
axis 'equal';  axis 'tight';     
colormap jet;
axis([xmin xmax ymin ymax]);
