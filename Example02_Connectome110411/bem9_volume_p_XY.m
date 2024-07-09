%   This script accurately computes and displays electric potential sampled
%   on a cross-section (transverse plane) via the FMM method with accurate
%   neighbor integration
%
%   Copyright SNM/WAW 2017-2024

%%  Load/prepare data
planeABCD = [0 0 1 -1e-3*Zm]; %  Analytical equation of the observation point plane for neighbor search acceleration

%% Define observation points in the cross-section
Ms = 300;
x = linspace(xmin, xmax, Ms);
y = linspace(ymin, ymax, Ms);
[X0, Y0]  = meshgrid(x, y);
clear pointsXY;
pointsXY(:, 1) = reshape(X0, 1, Ms^2);
pointsXY(:, 2) = reshape(Y0, 1, Ms^2);  
pointsXY(:, 3) = Z*ones(1, Ms^2);

%% Find the potential at each observation point
tic
pointsXY       = pointsXY;     % Convert back to m
Ppri           = zeros(Ms*Ms, 1);
Psec           = zeros(Ms*Ms, 1);
[~, Ppri]      = bemf3_inc_field_electric_plain_dipoles(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, 1e-3*pointsXY, 1e-2*prec);
R = 2;  %   precise integration
Psec           = bemf5_volume_field_potential(1e-3*pointsXY, c, P, t, Center, Area, normals, R, planeABCD, prec);
Ptotal         = Ppri + Psec;   
disp([newline 'Potential calculated in ' num2str(toc) ' s']);

%% Plot the potential in the cross-section
figure;
% Potential contour plot
temp      = Ptotal;
th1 =  +1e-4*max(+temp);             %   in V
th2 =  -1e-4*max(-temp);             %   in V
levels      = 20;
bemf2_graphics_vol_field_log(temp, th1, th2, levels, x, y);
xlabel('Distance x, mm');
ylabel('Distance y, mm');
title(strcat('Potential (V), ', label, '-in the transverse plane'));

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
grid on; set(gcf,'Color','White');