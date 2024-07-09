%   Prepare cross-section observation planes for subsequent visualizations
%   Plot contours of intersected compartments

%% Parameters
Xm      = X;        % in millimeters
Ym      = Y;        % in millimeters
Zm      = Z;        % in millimeters
delta   = 45;       % in millimeters
xmin = Xm-delta;
xmax = Xm+delta;
ymin = Ym-delta;
ymax = Ym+delta;
zmin = Zm-delta;
zmax = Zm+delta;

%%  Setup compartmental colors
%   SKIN
i = 1; 
color(i, :) = [0 0 0]; 
%   BONE
i = 2; 
color(i, :) = [0.6 0.6 0.6]; 
%   CSF
i = 3; 
color(i, :) = [1 0.5 0.0]; 
%   GM
i = 4;
color(i, :) = [1 0.75 0.65];
%   WM
i = 5;
color(i, :) = [0 0.7 0.7];
%   VENTRICLES
i = 6; 
color(i, :) = [1 0.75 0.65]; 
%   EYES
i = 7; 
color(i, :) = [1 0.75 0.65]; 

%% Compute contours of intersected tissues in the XY cross-section
tissues = length(tissue);
PofXY = cell(tissues, 1);   %   intersection nodes for a tissue
EofXY = cell(tissues, 1);   %   edges formed by intersection nodes for a tissue
TofXY = cell(tissues, 1);   %   intersected triangles
NofXY = cell(tissues, 1);   %   normal vectors of intersected triangles
countXY = [];               %   number of every tissue present in the slice
for m = 1:tissues 
    [Pi, ti, polymask, flag] = meshplaneintXY(PS{m}, tS{m}, eS{m}, TriPS{m}, TriMS{m}, Zm);
    if flag % intersection found                
        countXY               = [countXY m];
        PofXY{m}            = Pi;               %   intersection nodes
        EofXY{m}            = polymask;         %   edges formed by intersection nodes
        TofXY{m}            = ti;               %   intersected triangles
        NofXY{m}            = nS{m}(ti, :);     %   normal vectors of intersected triangles        
    end
end

%% Compute contours of intersected tissues in the XZ cross-section
tissues = length(tissue);
PofXZ = cell(tissues, 1);   %   intersection nodes for a tissue
EofXZ = cell(tissues, 1);   %   edges formed by intersection nodes for a tissue
TofXZ = cell(tissues, 1);   %   intersected triangles
NofXZ = cell(tissues, 1);   %   normal vectors of intersected triangles
countXZ = [];               %   number of every tissue present in the slice
for m = 1:tissues 
    [Pi, ti, polymask, flag] = meshplaneintXZ(PS{m}, tS{m}, eS{m}, TriPS{m}, TriMS{m}, Ym);
    if flag % intersection found                
        countXZ               = [countXZ m];
        PofXZ{m}            = Pi;               %   intersection nodes
        EofXZ{m}            = polymask;         %   edges formed by intersection nodes
        TofXZ{m}            = ti;               %   intersected triangles
        NofXZ{m}            = nS{m}(ti, :);     %   normal vectors of intersected triangles        
    end
end

%% Compute contours of intersected tissues in the YZ cross-section
tissues = length(tissue);
PofYZ = cell(tissues, 1);   %   intersection nodes for a tissue
EofYZ = cell(tissues, 1);   %   edges formed by intersection nodes for a tissue
TofYZ = cell(tissues, 1);   %   intersected triangles
NofYZ = cell(tissues, 1);   %   normal vectors of intersected triangles
countYZ = [];               %   number of every tissue present in the slice
for m = 1:tissues 
    [Pi, ti, polymask, flag] = meshplaneintYZ(PS{m}, tS{m}, eS{m}, TriPS{m}, TriMS{m}, Xm);
    if flag % intersection found                
        countYZ               = [countYZ m];
        PofYZ{m}            = Pi;               %   intersection nodes
        EofYZ{m}            = polymask;         %   edges formed by intersection nodes
        TofYZ{m}            = ti;               %   intersected triangles
        NofYZ{m}            = nS{m}(ti, :);     %   normal vectors of intersected triangles        
    end
end

%% Plot tissue contours in the XY cross-section
figure;

% Display the contours
for m = countXY
    edges           = EofXY{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofXY{m}(:, 1);       %   this is for the contour  
    points(:, 2)    = +PofXY{m}(:, 2);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
patch([xmin xmin xmax xmax],[ymin ymax ymax ymin], 'c', 'FaceAlpha', 0.35);
% General settings
title(['Cross-section in the transverse plane at Zm = ' num2str(Zm) ' m']);
xlabel('x, m'); ylabel('y, m');
axis 'equal';  axis 'tight'; 
set(gcf,'Color','White');
hold on;
bemf1_graphics_dipole(1e3*strdipolePplus, 1e3*strdipolePminus, 1e3*strdipoleCurrent, 1);

%% Plot tissue contours in the XZ cross-section
figure;

% Display the contours
for m = countXZ
    edges           = EofXZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofXZ{m}(:, 1);       %   this is for the contour  
    points(:, 2)    = +PofXZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
patch([xmin xmin xmax xmax], [zmin zmax zmax zmin], 'c', 'FaceAlpha', 0.35);
% General settings
title(['Cross-section in the coronal plane at Ym = ' num2str(Ym) ' mm']);
xlabel('x, mm'); ylabel('z, mm');
axis 'equal';  axis 'tight'; 
set(gcf,'Color','White');
hold on;
bemf1_graphics_dipole(1e3*strdipolePplus, 1e3*strdipolePminus, strdipoleCurrent, 2);

%% Plot tissue contours in the YZ cross-section
figure;

% Display the contours
for m = countYZ
    edges           = EofYZ{m};              %   this is for the contour
    points          = [];
    points(:, 1)    = +PofYZ{m}(:, 2);       %   this is for the contour  
    points(:, 2)    = +PofYZ{m}(:, 3);       %   this is for the contour
    patch('Faces', edges, 'Vertices', points, 'EdgeColor', color(m, :), 'LineWidth', 2.0);    %   this is contour plot
end
patch([ymin ymin ymax ymax], [zmin zmax zmax zmin], 'c', 'FaceAlpha', 0.35);
% General settings
title(['Cross-section in the sagittal plane at Xm = ' num2str(Xm) ' mm']);
xlabel('y, mm'); ylabel('z, mm');
axis 'equal';  axis 'tight'; 
set(gcf,'Color','White');
hold on;
bemf1_graphics_dipole(1e3*strdipolePplus, 1e3*strdipolePminus, strdipoleCurrent, 3);
