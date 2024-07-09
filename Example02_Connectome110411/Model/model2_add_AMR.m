%   This is a mesh processor script: it adds adaptive b-refinement to the
%   model mesh based around the given target point (e.g., dipole position)
%
%   Copyright SNM/WAW 2017-2024

%%  Refinement parameters
method      = 'manifold';       %   default
TAUBIN      = 'yes';            %   Taubin low-pass filtering: yes/no
AMRSTEPS    = 4;                %   Number of AMR steps
factor      = 8;                %   Parameter of the refinement cost function
prec        = 1e-02;            %   FMM precision (for first approximation of IE)

%%  Initialize first approximation of integral equation and compute MeanAbsC
EpriC       = cell(Compartments, 1);
bC          = cell(Compartments, 1);
CC          = cell(Compartments, 1);

tic
[Epri, ~]         = bemf3_inc_field_electric_plain_dipoles(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, Center, prec);
for m = 1:length(name)
    EpriC{m}        = Epri(Indicator==m, :);
    bC{m}           = 2*(Contrast(m)*sum(normalsC{m}.*EpriC{m}, 2));    %   Right-hand side of the BEM-FMM equation   
    CC{m}           = bC{m}.*AreaC{m};                                  %   First approximation for total charge per facet
end
MeanAbsC = mean(abs(cell2mat(CC))); 
PrimaryFieldTime = toc

%% Compartment by compartment refinement loop
tic
CostFunction = zeros(length(tissue), AMRSTEPS);
parfor m = 1:length(tissue)
    if m == 1
        continue;                                           %  do not refine skin
    end
    for n = 1:AMRSTEPS
        refine      = find(abs(CC{m})>factor*MeanAbsC);          %  indexes into faces to be refined
        CostFunction(m, n)   = std(CC{m})/(factor*MeanAbsC);     %   rows are the iterations
        if ~isempty(refine)        
            [ref, rem, indexref] = meshrefmanifold(PC{m}, tC{m}, normalsC{m}, refine);  
            if ~isempty(ref)
                if strcmp(TAUBIN, 'yes')
                    ref             = meshtaubin(ref);                  %   apply Taubin smoothing
                end
                tC{m}           = [rem.t; ref.t+size(rem.P, 1)];        %   unite both meshes into structure ref
                PC{m}           = [rem.P; ref.P];
                normalsC{m}     = [rem.normals; ref.normals];
                [PC{m}, tC{m}]  = fixmesh(PC{m}, tC{m});
                AreaC{m}        = meshareas(PC{m}, tC{m});
                %   Update total charges by using old results
                EpriC{m}(indexref, :) = [];                             %   delete results from refined facets
                Centertemp   = meshtricenter(ref.P, ref.t);
                [Etemp, ~]   = bemf3_inc_field_electric_plain_dipoles(strdipolePplus, strdipolePminus, strdipolesig, strdipoleCurrent, Centertemp, prec);
                EpriC{m}     = [EpriC{m}; Etemp];                 
                bC{m}        = 2*(Contrast(m)*sum(normalsC{m}.*EpriC{m}, 2));    %   Right-hand side of the BEM-FMM equation   
                CC{m}        = bC{m}.*AreaC{m}; 
            end
        end
    end
    [tC{m}]  = meshreorient(PC{m}, tC{m}, normalsC{m});
end  
CostFunction
bRefinementTime = toc

%% Reconstruct the combined locally refined model
tic
for m = 1:length(name)
    temp            = ones(size(tC{m}, 1), 1);
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
end
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
CombinedModelTime = toc

%%  Fix triangle orientation (just in case, optional)
t = meshreorient(P, t, normals);



  


