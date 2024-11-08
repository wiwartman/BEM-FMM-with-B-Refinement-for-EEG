clear all; %#ok
% First run: solution with first selection of AMR
ext_AMRSTEPS = 4; %#ok
ext_k = 6; %#ok
wrapper_execute;
soln_amr1.SurfaceP = PtotSkin;
soln_amr1.SurfaceB = Btotal;

% Second run: solution with second selection of AMR
ext_AMRSTEPS = 6; 
ext_k = 6;
wrapper_execute;
soln_amr2.SurfaceP = PtotSkin;
soln_amr2.SurfaceB = Btotal;

% Calculate and display 2-norm error between solutions
diffP = norm(soln_amr1.SurfaceP - soln_amr2.SurfaceP)/norm(soln_amr2.SurfaceP);
diffB = norm(soln_amr1.SurfaceB - soln_amr2.SurfaceB)/norm(soln_amr2.SurfaceB);

disp(['DiffP: ' num2str(diffP*100) '%']);
disp(['DiffB: ' num2str(diffB*100) '%']);