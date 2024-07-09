%   This script accurately computes and displays magnetic fields sampled on
%   a cross-section (transverse plane) via the FMM method with accurate
%   neighbor integration
%
%   Copyright SNM/WAW 2018-2024

%%  Load/prepare data
planeABCD = [0 0 1 -1e-3*Zm]; %  Analytical equation of the observation point plane for neighbor search acceleration

%%  Post processing parameters
component   = 4;        %   field component to be plotted (1, 2, 3 or x, y, z, or 4 - total) 
temp        = ['x' 'y' 'z' 't'];
label       = temp(component);

%%  Define observation points in the cross-section, mm
Ms = 200;
x = linspace(xmin, xmax, Ms);
y = linspace(ymin, ymax, Ms);
[X0, Y0]  = meshgrid(x, y);
clear pointsXY;
pointsXY(:, 1) = reshape(X0, 1, Ms^2);
pointsXY(:, 2) = reshape(Y0, 1, Ms^2);  
pointsXY(:, 3) = Zm*ones(1, Ms^2);

%% Find the B-field at each observation point
tic
pointsXY       = pointsXY;     % Convert back to m
Bpri           = zeros(Ms*Ms, 3);
Bsec           = zeros(Ms*Ms, 3);
difference     = condin - condout;
Bpri           = bemf3_inc_field_magnetic(strdipolemvector, strdipolemcenter, strdipolemstrength, 1e-3*pointsXY, mu0, 1e-2*prec);                                     
R = 2;         %   precise integration            
Bsec           = bemf5_volume_field_magnetic(1e-3*pointsXY, Ptot, P, t, Center, Area, normals, difference, mu0, R, 1, planeABCD, prec);
Btotal         = Bpri + Bsec;     
disp([newline 'B-field calculated in ' num2str(toc) ' s']);

%%  Plot the B-field in the cross-section
figure;
% B-field contour plot
if component == 4
    temp      = 1e12*abs(sqrt(dot(Btotal, Btotal, 2)));
else
    temp      = 1e12*abs(Btotal(:, component));
end
th1 = +1e-3*max(temp);              %   in T
th2 = 1*min(temp);                  %   in T
levels      = 20;
bemf2_graphics_vol_field_log(temp, th1, th2, levels, x, y);
xlabel('Distance x, mm');
ylabel('Distance y, mm');
title(strcat('B-field pT, ', label, '-component in the transverse plane.'));

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

%   General settings 
axis 'equal';  axis 'tight';     
colormap jet; 
axis([xmin xmax ymin ymax]);
grid on; set(gcf,'Color','White');