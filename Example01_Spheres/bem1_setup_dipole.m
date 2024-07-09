%   This script creates the base dipole (single)
%
%   Copyright SNM/WAW 2018-2020

%%   Single dipole example
I0 = 1e-6;              %   source current, A

%  Horizontal short dipole (0.04 mm length)
strdipolePplus      = [+0.00002 0.000 0.0760];  %   in m
strdipolePminus     = [-0.00002 0.000 0.0760];  %   in m

%   Rotate dipole if necessary
theta = 0.0;
strdipolePplus  = meshrotate2(strdipolePplus, [1 0 0], theta);
strdipolePminus = meshrotate2(strdipolePminus, [1 0 0], theta);
dlength         = norm(strdipolePplus - strdipolePminus)

normalization          = 1e-7/(I0*dlength); %    Dipole moment of 0.1 uA*m
strdipolesig           = [cond(4) cond(4)]';    % Change here
strdipoleCurrent       = normalization*[+I0 -I0]';
M                      = size(strdipolePplus, 1);

%%   Magnetic dipole subdivision (optional)
D = 1;                        %   number of smaller subdipoles
strdipolemvector   = zeros(D*M, 3);
strdipolemcenter   = zeros(D*M, 3);
strdipolemstrength = zeros(D*M, 1);
for m = 1:M
    temp = (1/D)*(strdipolePplus(m, :) - strdipolePminus(m, :));
    for d = 1:D 
        arg = d+D*(m-1);
        strdipolemvector(arg, :)     = temp;
        strdipolemcenter(arg, :)     = strdipolePminus(m, :) + (d-1/2)*temp;
        strdipolemstrength(arg, :)   = strdipoleCurrent(m);                  
    end
end

%%  Plot and check correct position   
R           = 0.015;
GM.t        = t(Indicator==4, :);
GM.P        = P;                 
[GM.P, GM.t] = fixmesh(GM.P, GM.t); 
GM.Center    = meshtricenter(GM.P, GM.t);

Ctr    = mean((strdipolePplus + strdipolePminus)/2, 1);

indexg1 = find( (GM.P(GM.t(:, 1), 1)-Ctr(1)).^2 + (GM.P(GM.t(:, 1), 2)-Ctr(2)).^2 + (GM.P(GM.t(:, 1), 3)-Ctr(3)).^2 < R^2);
indexg2 = find( (GM.P(GM.t(:, 2), 1)-Ctr(1)).^2 + (GM.P(GM.t(:, 2), 2)-Ctr(2)).^2 + (GM.P(GM.t(:, 2), 3)-Ctr(3)).^2 < R^2);
indexg3 = find( (GM.P(GM.t(:, 3), 1)-Ctr(1)).^2 + (GM.P(GM.t(:, 3), 2)-Ctr(2)).^2 + (GM.P(GM.t(:, 3), 3)-Ctr(3)).^2 < R^2);
indexg  = intersect(intersect(indexg1, indexg2), indexg3); 

% % Plot dipole below GM
f1 = figure;
str.EdgeColor = 'k'; str.FaceColor = 'c'; str.FaceAlpha = 1.0; 
bemf2_graphics_base(GM.P, GM.t(indexg, :), str);
bemf1_graphics_dipole(strdipolePplus, strdipolePminus, strdipoleCurrent, 4) 
axis 'equal';  axis 'tight';   
daspect([1 1 1]);
set(gcf,'Color','White');
xlabel('x'); ylabel('y'); zlabel('z');
view(0, 0); 
