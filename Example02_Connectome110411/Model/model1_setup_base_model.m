%   This is a mesh processor script: it loads model and the dipole
%   parameters into workspace and performs basic preparatory operations
%
%   Copyright SNM/WAW 2017-2024

%%  Setup path to engine
if ~isunix
    s = pwd; addpath(strcat(s(1:end-6), '\Engine'));
else
    s = pwd; addpath(strcat(s(1:end-6), '/Engine'));
end

%%  Setup parallel port
numThreads      = 12;      %   number of cores to be used
tic
ppool           = gcp('nocreate');
if isempty(ppool)
    tic
    parpool(numThreads);
    disp([newline 'Started parallel pool in ' num2str(toc) ' s']);
end

%%  Load dipole data and define local and global solutions
Radius          = 0.04;         %   Sphere radius for a local/two level solution, m
load('110411_dip1_data.mat');
X = 1e3*0.5*(strdipolePplus(1) + strdipolePminus(1));   % mm
Y = 1e3*0.5*(strdipolePplus(2) + strdipolePminus(2));   % mm
Z = 1e3*0.5*(strdipolePplus(3) + strdipolePminus(3));   % mm
R = 0.01;

%% Load compartment names and assign conductivities
index_name  = 'tissue_index.txt';
condambient = 0;
[name, tissue, cond, enclosingTissueIdx] = tissue_index_read(index_name);
[Contrast, Condin, Condout]              = assign_initial_conductivities(cond, condambient, enclosingTissueIdx);

%%  Construct compartmental meshes and their parameters
Compartments    = length(name);
PC              = cell(Compartments, 1);    %   for individual compartments
tC              = cell(Compartments, 1);    %   for individual compartments
ttemp           = cell(Compartments, 1);    %   for individual compartments
normalsC        = cell(Compartments, 1);    %   for individual compartments
CenterC         = cell(Compartments, 1);    %   for individual compartments
AreaC           = cell(Compartments, 1);    %   for individual compartments
contrastC       = cell(Compartments, 1);    %   for individual compartments
condinC         = cell(Compartments, 1);    %   for individual compartments
condoutC        = cell(Compartments, 1);    %   for individual compartments
IndicatorC      = cell(Compartments, 1);    %   for individual compartments 
tS              = cell(Compartments, 1);    %   for cross-sections only
nS              = cell(Compartments, 1);    %   for cross-sections only
eS              = cell(Compartments, 1);    %   for cross-sections only
TriPS           = cell(Compartments, 1);    %   for cross-sections only
TriMS           = cell(Compartments, 1);    %   for cross-sections only
PS              = cell(Compartments, 1);    %   for cross-sections only

tic
parfor m = 1:length(name)
    TR              = stlread(name{m});
    PC{m}           = 1e-3*TR.Points; %  only if the original data were in mm!
    tC{m}           = TR.ConnectivityList;
    temp            = ones(size(tC{m}, 1), 1);
    normalsC{m}     = meshnormals(PC{m}, tC{m}); 
    tC{m}           = meshreorient(PC{m}, tC{m}, normalsC{m});
    CenterC{m}      = meshtricenter(PC{m}, tC{m});
    AreaC{m}        = meshareas(PC{m}, tC{m});  
    contrastC{m}    = Contrast(m)*temp;
    condinC{m}      = Condin(m)*temp;
    condoutC{m}     = Condout(m)*temp;
    IndicatorC{m}   = m*temp;
    PS{m}           = 1e3*PC{m}; 
    tS{m}           = tC{m};
    nS{m}           = normalsC{m}; 
    [eS{m}, TriPS{m}, TriMS{m}] = mt(tS{m});
    disp(['Successfully loaded file [' name{m} ']']);
end
ModelLoadTime = toc

%%  Construct combined mesh for all compartments
tic
P           = cell2mat(PC);
shift       = 0;
for m = 1:length(name)
    ttemp{m} = tC{m} + shift;
    shift = shift + size(PC{m}, 1);
end
t           = cell2mat(ttemp);
normals     = cell2mat(normalsC);
Area        = cell2mat(AreaC);
Center      = cell2mat(CenterC);
contrast    = cell2mat(contrastC);
condin      = cell2mat(condinC);
condout     = cell2mat(condoutC);
Indicator   = cell2mat(IndicatorC);

%%  Create local mesh for a local AMR
D           = Center - repmat(strdipolemcenter, size(Center, 1), 1);
DIST        = sqrt(sum(D.*D, 2));
indexmain   = DIST<Radius;
L.t         = t(indexmain, :);
L.normals   = normals(indexmain, :);
L.Center    = Center(indexmain, :);
L.Area      = Area(indexmain);
L.Indicator = Indicator(indexmain);
L.condin    = condin(indexmain);
L.condout   = condout(indexmain);
L.contrast  = contrast(indexmain);
CombinedModelTime = toc