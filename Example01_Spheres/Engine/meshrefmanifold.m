function [ref, rem] = meshrefmanifold(obj, refine)
%   This script performs manifold mesh refinement using
%   barycentric triangle subdivision 1:4 and border subdivision 1:2
%   INPUTS
%   obj.P       - contains nodes of the mesh to be refined
%   obj.t       - set of facets of the mesh to be refined
%   obj.normals - normal vectors of the mesh to be refined
%   refine      - indexes of facets to be refined
%   OUTPUTS
%   ref.P       - nodes of the refined piece(s) of the mesh
%   ref.t       - facets of the refined piece(s) of the mesh
%   ref.normals - normal vectors of the refined piece(s) of the mesh
%   rem.P       - nodes of the remaining (non-refined) portion of the mesh
%   rem.P and ref.P intersect (share the same border nodes)
%   rem.t       - facets of the remaining (non-refined) portion of the mesh
%   rem.normals - normal vectors of the remaining (non-refined) portion of the mesh
%   Both ref and rem can be united.
%   SNM 2024

    %%  Return if a small dataset is used
    if length(refine)<10
        ref = [];
        rem = obj; 
        %   disp('no refinement is made')
        return;
    end

    %%   Find topological neighbors for the entire mesh
    DT = triangulation(obj.t, obj.P); 
    tneighbor = neighbors(DT);

    %% Expand dataset by including all non-refined facets with 2 or 3 neighbors to be refined
    %   After this block, any facet to be refined may contact no more than
    %   one face to be non-refined
    k = 1; warning off;
    while 1
        %  For every face to be non-refined find if the number of neighbors to be refined is >=2
        nonrefined  = setdiff([1:size(obj.t, 1)]', refine);
        DTN         = triangulation(obj.t(nonrefined, :), obj.P); 
        tneighborN  = neighbors(DTN);
        temp1       = isnan(tneighborN(:, 1));
        temp2       = isnan(tneighborN(:, 2));
        temp3       = isnan(tneighborN(:, 3));
        temp        = (temp1 + temp2 + temp3) >=2;
        if sum(temp) == 0
            break;
        end
        refine      = unique([refine; nonrefined(temp)], 'stable');
        k = k + 1;
    end
    warning on;

    %% Reduce dataset by excluding all facets to be refined with only 0 or 1 neighbors to be refined
    k = 1; warning off;
    while 1
        %   Return for a dataset bein too small
        if length(refine)<10
            ref = [];
            rem = obj; 
            %   disp('no refinement is made')
            return;
        end
        DTR         = triangulation(obj.t(refine, :), obj.P); 
        tneighborR  = neighbors(DTR);
        temp1       = isnan(tneighborR(:, 1));
        temp2       = isnan(tneighborR(:, 2));
        temp3       = isnan(tneighborR(:, 3));
        temp        = (temp1 + temp2 + temp3) >=2;
        indicator   = (~temp1 + ~temp2 + ~temp3)';
        if sum(temp) == 0
            break;
        end
        refine(temp)  = [];
        k = k + 1;
    end
    warning on;

    if length(refine)<10
        ref = [];
        rem = obj; 
        %   disp('no refinement is made')
        return;
    end

    %%  Sort all facets to be refined with regard to their neighbors to be refined: either 3 or 2
    index2      = find(indicator==2);
    index3      = find(indicator==3);
    refine      = refine([index3 index2]);
    indicator   = indicator([index3 index2]);
    F2          = length(index2);
    F3          = length(index3);

    %%  Define a neigbor sticking out for every border face (it will be divided into two facets)
    outerneigh = zeros(F2, 1);
    for m = 1:F2
        temp = setdiff(tneighbor(refine(F3+m), :), refine);
        outerneigh(m) = temp; % must always exist
    end

    %%  Construct a star set by including all outer neighbors (all local)
    P           = obj.P;
    t           = [obj.t(refine, :); obj.t(outerneigh, :)];
    t           = sort(t, 2);
    normals     = [obj.normals(refine, :); obj.normals(outerneigh, :)];
    R0          = length(refine);
    F0          = length(outerneigh);

    %%  Find all edges of the mesh to be refined excluding border edges
    edgesdup   = [t(:,[1, 2]); t(:,[1, 3]); t(:,[2, 3])];         %  All edges duplicated
    e1         = sort(edgesdup, 2);                               %  All edges now in same order 
    [edges,eia(:,1)] = unique(e1, 'rows');                        %  eia(:,1) contains the row in e1 that was chosen for edges
    [~,eia(:,2)] = unique(e1, 'last', 'rows');                    %  eia(:,2) contains the other row in e1 for the same edges
    e2t        = mod(eia-1,size(t,1))+1;                          %  since "edgesdup" contains three copies of t, mod makes the translation.
    e2t        = sort(e2t,2);                                     %  Put the lowest index triangle into TriP.
    e2tnondup = e2t(:,1) ~= e2t(:,2);                             %  Edges that are attached to only one triangle should be removed
    e2t = e2t(e2tnondup,:);
    edges = edges(e2tnondup,:);
    TriP       = e2t(:,1);
    TriM       = e2t(:,2);
    EdgesTotal = size(edges, 1);
    edgeSum    = sum(edges, 2);
    VerP       = sum(t(TriP, :), 2) - edgeSum; % Sum is much faster than setdiff in loop
    VerM       = sum(t(TriM, :), 2) - edgeSum;

    %%  Eliminate all edges between neighboring faces 
    nokeep  = (TriP>R0)&(TriM>R0);
    edges(nokeep, :)    = [];
    TriP(nokeep, :)     = [];
    TriM(nokeep, :)     = [];
    VerP(nokeep, :)     = [];
    VerM(nokeep, :)     = [];
    EdgesTotal          = size(edges, 1);
    
    %%  Plain subdivision of edges
    ecenter     = 1/2*(P(edges(:, 1), :) + P(edges(:, 2), :));
    newpoints   = ecenter;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Loop over facets F3 and F2 - find thee edge indexes (1 2/1 3/2 3)
    EdgeIndexes = zeros(F3+F2, 3);
    for m = 1:F3+F2
        temp = sum(abs(edges - repmat(t(m, [1 2]), EdgesTotal, 1)), 2);
        EdgeIndexes(m, 1) = find(temp==0);
        temp = sum(abs(edges - repmat(t(m, [1 3]), EdgesTotal, 1)), 2);
        EdgeIndexes(m, 2) = find(temp==0);
        temp = sum(abs(edges - repmat(t(m, [2 3]), EdgesTotal, 1)), 2);
        EdgeIndexes(m, 3) = find(temp==0);
    end

    %%  Loop over facets F3 and F2 - construct points to add
    P12 = zeros(F3+F2, 3);
    P13 = zeros(F3+F2, 3);
    P23 = zeros(F3+F2, 3);
    for m = 1:F3+F2
        P12(m, :) = newpoints(EdgeIndexes(m, 1), :);
        P13(m, :) = newpoints(EdgeIndexes(m, 2), :);
        P23(m, :) = newpoints(EdgeIndexes(m, 3), :);
    end

    %%  Define initial face set F3 and F2
    t32         = t(1:F3+F2, :);
    normals32   = normals(1:F3+F2, :);

    %% Introduce a duplicated set of nodes for set t32 to be refined 
    Pnew = [P; P12; P13; P23]; % size(P, 1), size(t32, 1), size(t32, 1), size(t32, 1)

    %%  Refine set t32
    indexp      = size(P, 1);  
    indext      = size(t32, 1);   
    array       = [1:indext]';
    tA(:, 1)    = array + indexp;                       %   12
    tA(:, 2)    = array + indexp + indext;              %   13    
    tA(:, 3)    = array + indexp + indext + indext;     %   23      
    tB(:, 1)    = t32(:, 1);                            %   1
    tB(:, 2)    = array + indexp;                       %   12    
    tB(:, 3)    = array + indexp + indext;              %   13 
    tC(:, 1)    = t32(:, 2);                            %   2
    tC(:, 2)    = array + indexp;                       %   12    
    tC(:, 3)    = array + indexp + indext + indext;     %   23 
    tD(:, 1)    = t32(:, 3);                            %   3
    tD(:, 2)    = array + indexp + indext;              %   13    
    tD(:, 3)    = array + indexp + indext + indext;     %   23        
    tnew = [tA; tB; tC; tD];
    %  Remove duplicated nodes
    [P32, t32]  = fixmesh(Pnew, tnew);  
    %  Update normal vectors
    basenormals = meshnormals(P32, t32);
    unitnormals = [normals32; normals32; normals32; normals32];
    SIGN        = sign(sum(basenormals.*unitnormals, 2));
    normals32   = repmat(SIGN, 1, 3).*basenormals;
    %   Do the rest
    t32         = meshreorient(P32, t32, normals32);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  Loop over border facets F0 - find edge index (one basis function)
    EdgeIndexes = zeros(F0, 4);
    for m = 1:F0
        temp1 = sum(abs(edges - repmat(t(m+R0, [1 2]), EdgesTotal, 1)), 2);
        temp1 = find(temp1==0);
        temp2 = sum(abs(edges - repmat(t(m+R0, [1 3]), EdgesTotal, 1)), 2);
        temp2 = find(temp2==0);
        temp3 = sum(abs(edges - repmat(t(m+R0, [2 3]), EdgesTotal, 1)), 2);
        temp3 = find(temp3==0);
        if ~isempty(temp1)
            EdgeIndexes(m, :) = [1 2 3 temp1];
        end
        if ~isempty(temp2)
            EdgeIndexes(m, :) = [1 3 2 temp2];
        end
        if ~isempty(temp3)
            EdgeIndexes(m, :) = [2 3 1 temp3];
        end
    end

    %%  Loop over border facets F0 - construct points to add
    P11 = zeros(F0, 3);
    for m = 1:F0
        P11(m, :) = newpoints(EdgeIndexes(m, 4), :);
    end

    %%  Loop over border facets F0 - construct facets to add
    Pnew        = [P; P11];    % size(P, 1), [F0, 3]
    PN          = size(P, 1);
    tnew        = zeros(2*F0, 3); %   subdivide into two
    normals10   = zeros(2*F0, 3); %   subdivide into two
    for m = 1:F0
        vertex1 = EdgeIndexes(m, 1);
        vertex2 = EdgeIndexes(m, 2);
        vertex3 = EdgeIndexes(m, 3);
        tnew(2*m-1, :) = [t(m+R0, vertex1) t(m+R0, vertex3) PN+m];
        tnew(2*m-0, :) = [t(m+R0, vertex2) t(m+R0, vertex3) PN+m];
        normals10(2*m-1, :) = normals(m+R0, :);
        normals10(2*m-0, :) = normals(m+R0, :);
    end

    %  Remove duplicated nodes   
    [P10, t10]  = fixmesh(Pnew, tnew);  
    %  Update normal vectors
    basenormals = meshnormals(P10, t10);
    unitnormals = normals10;
    SIGN        = sign(sum(basenormals.*unitnormals, 2));
    normals10   = repmat(SIGN, 1, 3).*basenormals;
    %   Do the rest
    t10         = meshreorient(P, t10, normals10);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%   Unite three sets of data
    %   SetNR (left nonrefined)
    P1 = obj.P;
    t1 = obj.t; 
    t1([refine; outerneigh], :) = []; 
    normals1 = obj.normals; 
    normals1([refine; outerneigh], :) = []; 

    %   SetR4 (refined 1:4, set32)
    P2 = P32;
    t2 = t32; 
    normals2 = normals32;  

    %   SetR2 (refined 1:2, set10)
    P3 = P10;
    t3 = t10; 
    normals3 = normals10;  

    ref.t           = [t1; t2+size(P1, 1); t3+size(P1, 1)+size(P2, 1)];
    ref.P           = [P1; P2; P3];
    ref.normals     = [normals1; normals2; normals3];
    [ref.P, ref.t]  = fixmesh(ref.P, ref.t);
    ref.t           = meshreorient(ref.P, ref.t, ref.normals);

    %%   Unite/prepare refined dataset first
    %   SetR4 (refined 1:4, set32)
    PA = P32;
    tA = t32; 
    normalsA = normals32;  
    %   SetR2 (refined 1:2, set10)
    PB = P10;
    tB = t10; 
    normalsB        = normals10;  
    tR              = [tA; tB+size(PA, 1)];
    PR              = [PA; PB];
    normalsR        = [normalsA; normalsB];
    [ref.P, ref.t]  = fixmesh(PR, tR);
    ref.normals     = normalsR;
    ref.t           = meshreorient(ref.P, ref.t, ref.normals);

    %%   Prepare remaining dataset next
    %   SetNR (left nonrefined)
    PN                                  = obj.P;
    tN                                  = obj.t; 
    tN([refine; outerneigh], :)         = []; 
    normalsN                            = obj.normals; 
    normalsN([refine; outerneigh], :)   = []; 
    [PN, tN]                            = fixmesh(PN, tN);
    tN                                  = meshreorient(PN, tN, normalsN);
    rem.P                               = PN;
    rem.t                               = tN;
    rem.normals                         = normalsN;
end