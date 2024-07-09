function [taubin]     = meshtaubin(obj)
%   Input to this function is a mesh to be smoothed
%   Input mesh does not have to be manifold
%   Border nodes (if any) will not be smoothed
%   INPUTS
%   obj.P       - contains nodes of the mesh to be smoothed
%   obj.t       - set of facets of the mesh to be smoothed
%   obj.normals - normal vectors of the mesh to be smoothed
%   OUTPUTS
%   taubin.P       - nodes of the smoothed mesh
%   taubin.t       - facets of the smoothed mesh
%   taubin.normals - normal vectors of the smoothed mesh
%   SNM 2024

    P       = obj.P;
    t       = obj.t;
    normals = obj.normals;
    
    %%  Find border nodes not to be moved
    TR          = triangulation(t, P);
    e           = edges(TR);
    [ft, fP]    = freeBoundary(TR);
    [C, ia]     = setdiff(P, fP, 'rows', 'stable');     %   indexes into non-border nodes
    ifP         = setdiff([1:size(P, 1)], ia);          %   indexes into border nodes (not to be moved)
    
    %%  Apply Taubin smoothing
    surfaceMeshIn = surfaceMesh(P, t);
    numIterations = 5;
    scaleFactor = [-0.62 0.6];
    surfaceMeshOut = smoothSurfaceMesh(surfaceMeshIn,numIterations,Method="Taubin",ScaleFactor=scaleFactor);
    taubin.P = surfaceMeshOut.Vertices;
    taubin.t = surfaceMeshOut.Faces;
    
    %%  Restore old border nodes
    taubin.P(ifP, :) = P(ifP, :); 
    
    %%   Compute normal vectors;
    taubin.normals      = meshnormals(taubin.P, taubin.t); 
    SIGN                = sign(sum(normals.*taubin.normals, 2));
    taubin.normals      = repmat(SIGN, 1, 3).*taubin.normals;
end
