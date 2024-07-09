%   This is a mesh processor script: it constructs a combined mesh of a
%   multi-object structure (for example, a head or a whole body)
%
%   Copyright SNM/WAW 2017-2024

clear all; %#ok<CLALL>

%%  Define EM constants
eps0        = 8.85418782e-012;  %   Dielectric permittivity of vacuum(~air) F/m
mu0         = 1.25663706e-006;  %   Magnetic permeability of vacuum(~air) H/m

%%  Add path (Windows/Linux commands are different)
s = pwd;
if(~isunix)
    slash = '\';
else
    slash = '/';
end
%   If a simulation has already been run with a different model,
%   that model's path may have been loaded, and the desired model's 
%   files would then be shadowed by the files on the old model path
warning off; rmpath(genpath(s)); warning on;
engine_path =   [s, slash, 'Engine']; 
addpath(engine_path);

%% Load tissue filenames and tissue display names from index file
index_name = 'tissue_index.txt';
[name, tissue, cond, enclosingTissueIdx] = tissue_index_read(index_name);

%%  Load tissue meshes and combine individual meshes into a single mesh
tic
PP = [];
tt = [];
nnormals = [];
Indicator = [];

%   Combine individual meshes into a single mesh
factor = [920 860 800 780];
for m = 1:length(name)
    load(name{m}); 
    P = P*factor(m);    %  for sphere of 0.1 in radius
    P = P*1e-3;         %  only if the original data were in mm!
    tt = [tt; t+size(PP, 1)];
    PP = [PP; P];
    nnormals = [nnormals; normals];    
    Indicator= [Indicator; repmat(m, size(t, 1), 1)];
    disp(['Successfully loaded file [' name{m} ']']);
end
t = tt;
P = PP;
normals = nnormals;
disp([newline 'Base meshes loaded in ' num2str(toc) ' s']);

%%  Fix triangle orientation (just in case, optional)
tic
t = meshreorient(P, t, normals);

%%   Process other mesh data
Center      = 1/3*(P(t(:, 1), :) + P(t(:, 2), :) + P(t(:, 3), :));  %   face centers
Area        = meshareas(P, t);  
disp([newline 'Triangle properties computed in ' num2str(toc) ' s']);

%%  Assign facet conductivity information
tic
condambient = 0.0; %   air
[contrast, condin, condout] = assign_initial_conductivities(cond, condambient, Indicator, enclosingTissueIdx);
disp([newline 'Initial conductivity contrasts assigned in ' num2str(toc) ' s']);

%%  Check for and process triangles that have coincident centroids
tic
disp('Checking combined mesh for duplicate facets ...');
[P, t, normals, Center, Area, Indicator, condin, condout, contrast] = ...
    clean_coincident_facets(P, t, normals, Center, Area, Indicator, condin, condout, contrast);
disp('Resolved all duplicate facets');
N           = size(t, 1);
disp([newline 'Duplicate facets resolved in ' num2str(toc) ' s']);
