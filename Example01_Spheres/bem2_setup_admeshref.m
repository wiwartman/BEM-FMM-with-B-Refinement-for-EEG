%   This script performs adaptive mesh refinement based on the first
%   appriximation for surface charge ( proprtional to b)
%
%   Copyright SNM/WAW 2021-2024

%%  Multiple refinement steps
method      = 'manifold';
TAUBIN      = 'yes';
AMRSTEPS    = 4; 


tic
for n = 1:AMRSTEPS
    %   Compute RHS and establish refinement criterion
    [Epri, ~] = bemf3_inc_field_electric(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, P, t, Center, Area, normals, 1e-2, 0, 0);
    b        = 2*(contrast.*sum(normals.*Epri, 2)); C = b.*Area; %  Right-hand side of the BEM-FMM equation    
    %   Assemble refinement
    PP          = [];
    tt          = [];
    nnormals    = [];
    nIndicator  = [];
    %% Compartment by compartment refinement loop
    for m = 1:length(tissue)
        index = Indicator==m; 
        % Refine selected compartment
        obj.P = P; obj.t = t(index, :); obj.normals = normals(index, :);
        [obj.P, obj.t] = fixmesh(P, obj.t);                           %   optional
        refine         = find(abs(C(Indicator==m))>mean(abs(C)));     %   indexes into faces to be refined
        CostFunction(m, n)   = std(C(Indicator==m))/mean(abs(C));     %   rows are the iterations  
        if strcmp(method, 'manifold')
            [ref, rem] = meshrefmanifold(obj, refine);
            %   unite both meshes into ref
            if ~isempty(ref)
                if strcmp(TAUBIN, 'yes')
                    ref             = meshtaubin(ref);                    %   apply Taubin smoothing
                end
                ref.t           = [rem.t; ref.t+size(rem.P, 1)];
                ref.P           = [rem.P; ref.P];
                ref.normals     = [rem.normals; ref.normals];
                [ref.P, ref.t]  = fixmesh(ref.P, ref.t);
            else
                ref = rem;
            end
        else
            ref = meshrefnonmanifold(obj, refine);
        end
        %   Output
        tt = [tt; ref.t+size(PP, 1)];
        PP = [PP; ref.P];
        nnormals = [nnormals; ref.normals];    
        nIndicator= [nIndicator; repmat(m, size(ref.t, 1), 1)];
    end  
    t = tt;
    P = PP;
    normals     = nnormals;
    Indicator   = nIndicator;
    Center      = 1/3*(P(t(:, 1), :) + P(t(:, 2), :) + P(t(:, 3), :));  %   face centers
    Area        = meshareas(P, t); 
    [contrast, condin, condout] = assign_initial_conductivities(cond, 0.0, Indicator, enclosingTissueIdx);
end
disp([newline 'Mesh refinement calculated in ' num2str(toc) ' s']);
CostFunction

%%  Fix triangle orientation (just in case, optional)
tic
t = meshreorient(P, t, normals);


%%  Check for and process triangles that have coincident centroids
tic
disp('Checking combined mesh for duplicate facets ...');
[P, t, normals, Center, Area, Indicator, condin, condout, contrast] = ...
    clean_coincident_facets(P, t, normals, Center, Area, Indicator, condin, condout, contrast);
disp('Resolved all duplicate facets');
N           = size(t, 1);
disp([newline 'Duplicate facets resolved in ' num2str(toc) ' s']);

%%   Find topological neighbors
tic
DT = triangulation(t, P); 
tneighbor = neighbors(DT);
% Fix cases where not all triangles have three neighbors
tneighbor = pad_neighbor_triangles(tneighbor);

%%  Graphics
for m = 1:length(tissue)
    figure;
    t0              = t(Indicator==m, :);   
    NumberOfTrianglesInShell = size(t0, 1);   
    patch('faces', t0, 'vertices', P, 'FaceColor', 'c', 'EdgeColor', 'k', 'FaceAlpha', 1.0);   
    axis equal; axis tight;
    set(gcf,'Color','White'); 
    axis(1e-3*[-95 95 -95 95 -95 95]);
    axis off;
    %light
end

  


