function ref = meshrefnonmanifold(obj, refine)
%   This script performs non-manifold mesh refinement using
%   barycentric triangle subdivision 1:4 only
%   INPUTS
%   obj.P       - contains nodes of the mesh to be refined
%   obj.t       - set of facets of the mesh to be refined
%   obj.normals - normal vectors of the mesh to be refined
%   refine      - indexes of facets to be refined
%   OUTPUTS
%   ref.P       - nodes of the refined mesh
%   ref.t       - facets of the refined mesh
%   ref.normals - normal vectors of the refined mesh
%   SNM 2021-2024

    %%  Return if a small dataset is used
    if length(refine)<10
        ref = obj;
        disp('no refinement is made')
        return;
    end

    P           = obj.P;
    t           = obj.t(refine, :);
    normals     = obj.normals(refine, :);

    %% Introduce new nodes
    %   New nodes 12 (between vertices 1 and 2) for every triangle
    P12 = 0.5*(P(t(:, 1), :) + P(t(:, 2), :));
    %   New nodes 13 (between vertices 1 and 3) for every triangle
    P13 = 0.5*(P(t(:, 1), :) + P(t(:, 3), :));
    %   New nodes 23 (between vertices 2 and 3) for every triangle
    P23 = 0.5*(P(t(:, 2), :) + P(t(:, 3), :));
    
    %% Introduce a duplicated set of nodes
    Pnew = [P; P12; P13; P23]; % size(P, 1), size(t, 1), size(t, 1), size(t, 1)
    
    %  Introduce the full set of triangles
    indexp      = size(P, 1);  
    indext      = size(t, 1);   
    array       = [1:indext]';
    
    tA(:, 1)    = array + indexp;                       %   12
    tA(:, 2)    = array + indexp + indext;              %   13    
    tA(:, 3)    = array + indexp + indext + indext;     %   23      
    
    tB(:, 1)    = t(:, 1);                              %   1
    tB(:, 2)    = array + indexp;                       %   12    
    tB(:, 3)    = array + indexp + indext;              %   13 
    
    tC(:, 1)    = t(:, 2);                              %   2
    tC(:, 2)    = array + indexp;                       %   12    
    tC(:, 3)    = array + indexp + indext + indext;     %   23 
    
    tD(:, 1)    = t(:, 3);                              %   3
    tD(:, 2)    = array + indexp + indext;              %   13    
    tD(:, 3)    = array + indexp + indext + indext;     %   23        
    
    tnew = [tA; tB; tC; tD];
    
    %  Remove duplicated nodes    
    [P, t]  = fixmesh(Pnew, tnew);  
    normals = [normals; normals; normals; normals]; 
    ref.P   = P;
    ref.t   = t;
    ref.normals = normals;

    obj.P                                   = [ref.P; obj.P];
    norefine                                = setdiff([1:length(obj.t)], refine);
    obj.t                                   = [ref.t; obj.t(norefine, :)+size(ref.P, 1)];
    obj.normals                             = [ref.normals; obj.normals(norefine, :)];     
    [obj.P, obj.t]                          = fixmesh(obj.P, obj.t);  
    ref                                     = obj; 
end