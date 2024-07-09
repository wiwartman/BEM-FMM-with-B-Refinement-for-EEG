%   This script computes the induced surface charge density for an
%   inhomogeneous multi-tissue object given the primary electric field, with
%   accurate neighbor integration
%
%   Copyright SNM/WAW 2017-2020

%%  Parameters of the iterative solution
iter         = 30;              %    Maximum possible number of iterations in the solution 
relres       = 1e-12;           %    Minimum acceptable relative residual 
weight       = 1/2;             %    Weight of the charge conservation law to be added (empirically found)

R       = 10;
flag    = 1;
[Epri, Ppri] = bemf3_inc_field_electric(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, P, t, Center, Area, normals, 1e-6, R, flag);
disp([newline 'Primary field calculated in ' num2str(toc) ' s']);
b        = 2*(contrast.*sum(normals.*Epri, 2));                         %  Right-hand side of the BEM-FMM equation

%%  GMRES iterative solution (native MATLAB GMRES is used)
h           = waitbar(0.5, 'Please wait - Running MATLAB GMRES'); 
tic
%   MATVEC is the user-defined function of c equal to the left-hand side of the matrix equation LHS(c) = b
MATVEC = @(c) bemf4_surface_field_lhs(c, Center, Area, contrast, normals, weight, EC);     
[c, its, resvec] = fgmres(MATVEC, b, relres, 'restart', iter, 'max_iters', 1, 'x0', b);
close(h);

figure
RESVEC = [];
for m = 1:size(resvec, 2)
    if m == size(resvec, 2)
        RESVEC = [RESVEC; resvec(1:its(2), m)];
    else
        RESVEC = [RESVEC; resvec(:, m)];
    end
end
semilogy(RESVEC, '-o'); grid on;
title('Relative residual of the iterative solution');
xlabel('Iteration number');
ylabel('Relative residual');

%%   Find surface electric potential
Padd = bemf4_surface_field_potential_accurate(c, Center, Area, PC);
Ptot = Ppri + Padd;     %   Continuous total electric potential at interfaces
