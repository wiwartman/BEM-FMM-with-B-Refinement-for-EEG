%   Parameters of BEM-FMM solution
%
%   Copyright SNM/WAW 2017-2024

eps0        = 8.85418782e-012;  %   Dielectric permittivity of vacuum(~air)
mu0         = 1.25663706e-006;  %   Magnetic permeability of vacuum(~air)
condambient = 0.0;              %   Air

%   Parameters for neighbor integrals
RnumberE        = 32;           %   Number of neighbor triangles for analytical integration of electric field - global solution
RnumberP        = 1;            %   Number of neighbor triangles for analytical integration of electric potential - global solution
L.RnumberE      = 32;           %   Number of neighbor triangles for analytical integration of electric field - local solution

%  Parameters of local iterative solution
Radius          = 0.04;         %   Sphere radius for a local GMRES solution
L.iter          = 30;           %    Maximum possible number of iterations in the solution 
L.relres        = 1e-6;         %    Minimum acceptable relative residual 
L.weight        = 1/2;          %    Weight of the charge conservation law to be added (empirically found)
%  Parameters of global iterative solution
iter            = 30;           %    Maximum possible number of iterations in the solution 
relres          = 5e-3;         %    Minimum acceptable relative residual 
weight          = 1/2;          %    Weight of the charge conservation law to be added (empirically found)

%   FMM Parameters
prec            = 1e-2;         %   Good for surfaces with neighbor integrals

%%  Create local mesh for a two-step GMRES
D      = Center - repmat(strdipolemcenter, size(Center, 1), 1);
DIST   = sqrt(sum(D.*D, 2));
index  = DIST<Radius;
L.t         = t(index, :);
L.normals   = normals(index, :);
L.Center    = Center(index, :);
L.Area      = Area(index);
L.Indicator = Indicator(index);
L.condin    = condin(index);
L.condout   = condout(index);
L.contrast  = contrast(index);






  


