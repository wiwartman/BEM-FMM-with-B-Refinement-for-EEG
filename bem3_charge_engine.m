%   This script computes the induced surface charge density for an
%   inhomogeneous multi-tissue object given the primary electric field, with
%   accurate neighbor integration
%
%   Copyright SNM/WAW 2017-2020


%%  Global accurate primary field
tic
gaussRadius = 2*R;
[Epri, Ppri] = bemf3_inc_field_electric_gauss_selective_dipoles(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, P, t, Center, Ctr, gaussRadius, prec);
PrimaryFieldTime = toc

%%  Local iterative solution
b        = 2*(L.contrast.*sum(L.normals.*Epri(index, :), 2));           %  Right-hand side of the BEM-FMM equation (local)
h        = waitbar(0.5, 'Please wait - Running MATLAB GMRES'); 
tic
%   MATVEC is the user-defined function of c equal to the left-hand side of the matrix equation LHS(c) = b
MATVEC = @(c) bemf4_surface_field_lhs(c, L.Center, L.Area, L.contrast, L.normals, L.weight, L.EC, prec);     
[c, ~, ~] = fgmres(MATVEC, b, L.relres, 'restart', L.iter, 'max_iters', 1, 'x0', b);
close(h);

%%  Append c (construct c0 for global solution)
C           = zeros(size(t, 1), 1);
C(index)    = c;
c0          = C; 

%%  Construct RHS for global solution
b        = 2*(contrast.*sum(normals.*Epri, 2));                         %  Right-hand side of the BEM-FMM equation
bprime   = b - bemf4_surface_field_lhs(c0, Center, Area, contrast, normals, weight, EC, prec); 

%%  Global iterative solution
h           = waitbar(0.5, 'Please wait - Running MATLAB GMRES'); 
tic
%   MATVEC is the user-defined function of cprime equal to the left-hand
%   side of the matrix equation LHS(cprime) = bprime
MATVEC = @(cprime) bemf4_surface_field_lhs(cprime, Center, Area, contrast, normals, weight, EC, prec);     
[cprime, its, resvec] = fgmres(MATVEC, bprime, relres, 'restart', iter, 'max_iters', 1, 'x0', bprime);
close(h);

%%  Add both solutions together
c = c0 + cprime;

%%   Find surface electric potential
tic
Padd = bemf4_surface_field_potential_accurate(c, Center, Area, PC, prec);
Ptot = Ppri + Padd;     %   Continuous total electric potential at interfaces
PotentialTime = toc
