function [problemDef,problemDef_cells,problemDef_limits,priors,controls] = RatParseClassToStructs_new(inputProblemDef,inputControls)

% Breaks up the classes into the relevant structures for inputting into C

% Put the extracted fields into a cell array...
% Structure of problemDef_cells array.
%
% {1} - inputProblemDef.contrastRepeatSLDs
%       {1 x nContrasts} array of cells
%       Each cell is {1 x 2 double}.
%
% {2} - inputProblemDef.allData
%       {1 x nContrasts} array of cells
%       Each cell is {Inf x 3 double}
%
% {3} - inputProblemDef.dataLimits
%       {1 x nContrasts} array of cells
%       Each cell is {1 x 2 double}
%
% {4} - inputProblemDef.simLimits
%       {1 x nContrasts} array of cells
%       Each cell is {1 x 2 double}
%
% {5} - inputProblemDef.contrastLayers
%       {1 x nContrasts} array of cells
%       Each cell is {1 x Inf double}
%
% {6} - inputProblemDef.layersDetails
%       {n x 1} array of cells
%       Each cell is (1 x 5 double}
%
% {7} - inputProblemDef.paramNames
%       {1 x nParams} array of cells
%       Each cell is {1 x Inf char}
%
% {8} - inputProblemDef.backgroundsNames
%       {1 x nBackgrounds} array of cells
%       Each cell is {1 x Inf char}
% 
% {9} - inputProblemDef.scalefactorNames
%       {1 x nScales} array of cells
%       Each cell is {1 x Inf char}
% 
% {10}- inputProblemDef.qzshiftNames
%       {1 x nShifts} array of cells
%       Each cell is {1 x Inf char}
% 
% {11}- inputProblemDef.nbairNames
%       {1 x nNba} array of cells
%       Each cell is {1 x Inf char}
% 
% {12}- inputProblemDef.nbsrNames
%       {1 x nNba} array of cells
%       Each cell is {1 x Inf char}
% 
% {13}- inputProblemDef.resolNames
%       {1 x nNba} array of cells
%       Each cell is {1 x Inf char}
 
% First parse the class to a structure variable.

inputStruct = inputProblemDef.toStruct();

%Start by removing the cell arrays
repeatLayers = inputStruct.contrastRepeatSLDs; %*****
allData = inputStruct.allData;
dataLimits = inputStruct.dataLimits;
simLimits = inputStruct.simLimits;
contrastLayers = inputStruct.contrastLayers;
layersDetails = inputStruct.layersDetails;
paramNames = inputStruct.paramNames;
paramPriors = inputStruct.paramPriors;
backsNames = inputStruct.backParNames;
backsPriors = inputStruct.backsPriors; % *****8 ToDo
sfNames = inputStruct.scalefactorNames;
shiftsNames = inputStruct.qzshiftNames; % TODO
shiftPriors = inputStruct.qzshiftPriors;
nbaNames = inputStruct.nbairNames;
nbaPriors = inputStruct.nbaPriors;
nbsNames = inputStruct.nbsubNames;
nbsPriors = inputStruct.nbsPriors;
resolNames = inputStruct.resolParNames;         % ******* ToDo
resolParPriors = inputStruct.resolParPriors;


% Pull out all the cell arrays (except priors) into one array
problemDef_cells{1} = repeatLayers;
problemDef_cells{2} = allData;
problemDef_cells{3} = dataLimits;
problemDef_cells{4} = simLimits;
problemDef_cells{5} = contrastLayers;

% Fix for cell array bug with custom layers
if strcmpi(inputStruct.modelType,'custom layers') || strcmpi(inputStruct.modelType,'Custom XY')
    for i = 1:length(problemDef_cells{5})
        problemDef_cells{5}{i} = 0;
    end
end

problemDef_cells{6} = layersDetails;
problemDef_cells{7} = paramNames;
problemDef_cells{8} = backsNames;             % ****** Todo
problemDef_cells{9} = sfNames;
problemDef_cells{10} = shiftsNames;
problemDef_cells{11} = nbaNames;
problemDef_cells{12} = nbsNames;
problemDef_cells{13} = resolNames;

% Put the priors into their own array
priors.paramPriors = paramPriors;
priors.backsPriors = backsPriors;
priors.resolPriors = resolParPriors;
priors.nbaPriors = nbaPriors;
priors.nbsPriors = nbsPriors;
priors.shiftPriors = shiftPriors;

%Split up the contrastBacks array
contrastBacks = inputStruct.contrastBacks;
for i = 1:length(contrastBacks)
    problemDef.contrastBacks(i) = contrastBacks{i}(1);
    problemDef.contrastBacksType(i) = contrastBacks{i}(2);
end
    
%Now make the limits array
for i = 1:length(inputStruct.paramConstr)
    problemDef_limits.params(i,:) = inputStruct.paramConstr{i};
end

for i = 1:length(inputStruct.backParconstr)
    problemDef_limits.backs(i,:) = inputStruct.backParconstr{i};
end

for i = 1:length(inputStruct.scalefactorConstr)
    problemDef_limits.scales(i,:) = inputStruct.scalefactorConstr{i};
end

for i = 1:length(inputStruct.qzshiftConstr)
    problemDef_limits.shifts(i,:) = inputStruct.qzshiftConstr{i};
end

for i = 1:length(inputStruct.nbairConstr)
    problemDef_limits.nba(i,:) = inputStruct.nbairConstr{i};
end

for i = 1:length(inputStruct.nbsubConstr)
    problemDef_limits.nbs(i,:) = inputStruct.nbsubConstr{i};
end

for i = 1:length(inputStruct.resolParConstr)
    problemDef_limits.res(i,:) = inputStruct.resolParConstr{i};
end

%Now remove all these fields from inputProblemDef
removedFields = {'contrastRepeatSLDs',...
    'allData',...
    'dataLimits',...
    'simLimits',...
    'contrastLayers',...
    'layersDetails',...
    'paramNames',...
    'backgroundNames',...
    'scalefactorNames',...
    'qzshiftNames',...
    'nbairNames',...
    'nbsubsNames',...
    'resolutionNames',...
    'paramConstr',...
    'backgroundConstr',...
    'scalefactorConstr',...
    'nbairConstr',...
    'nbsubConstr',...
    'resolutionConstr'};

%Make the problemDef structure from the bits left.....

% *************************************************************************
% NOTE - not using the more complicated background and resolution
% definitions for now - instead use the background names and backsPar
% values.... fix this next
% **********************************************************************8


problemDef.TF = 'standardTF'; %for now....
problemDef.resample = inputStruct.resample;
problemDef.dataPresent = inputStruct.dataPresent;
problemDef.numberOfContrasts = inputStruct.numberOfContrasts;
problemDef.geometry = inputStruct.geometry;
%problemDef.contrastBacks = contrastBacks;
problemDef.contrastShifts = inputStruct.contrastShifts;
problemDef.contrastScales = inputStruct.contrastScales;
problemDef.contrastNbas = inputStruct.contrastNbas;
problemDef.contrastNbss = inputStruct.contrastNbss;
problemDef.contrastRes = inputStruct.contrastRes;
problemDef.backs = inputStruct.backParVals; %inputStruct.backgrounds;             % **** note vackPar workaround (todo) ****
problemDef.shifts = inputStruct.qzshifts;
problemDef.sf = inputStruct.scalefactors;
problemDef.nba = inputStruct.nbairs;
problemDef.nbs = inputStruct.nbsubs;
problemDef.res = inputStruct.resolPars; %inputStruct.resolutions;           % **** note resolPar workaround (todo) ****          
problemDef.params = inputStruct.params;
problemDef.numberOfLayers = inputStruct.numberOfLayers;
problemDef.modelType = inputStruct.modelType;

if isfield(inputStruct,'modelFilename')
    if ~isempty(inputStruct.modelFilename)
        [path,fname,extension] = fileparts(inputStruct.modelFilename);
    else
        fname = '';
        path = pwd;
    end
else
    fname = '';
    path = pwd;
end
    
problemDef.modelFilename = fname;
problemDef.path = path;

if isfield(inputStruct,'modelLanguage')
    if ~isempty(inputStruct.modelLanguage)
        problemDef.modelLanguage = inputStruct.modelLanguage;
    else
        problemDef.modelLanguage = 'matlab';
    end
else
    problemDef.modelLanguage = '';
end
    
problemDef.fitpars = [];
problemDef.otherpars = [];
problemDef.fitconstr = [];
problemDef.otherconstr = [];

%Now deal with the controls class
controls.para = inputControls.parallel;
controls.proc = inputControls.procedure;
controls.display = inputControls.display;
controls.tolX = inputControls.tolX;
controls.tolFun = inputControls.tolFun;
controls.maxFunEvals = inputControls.maxFunEvals;
controls.maxIter = inputControls.maxIter;
controls.populationSize = inputControls.populationSize;
controls.F_weight = inputControls.F_weight;
controls.F_CR = inputControls.crossoverProbability;
controls.VTR = inputControls.targetValue;
controls.numGen = inputControls.numGenerations;
controls.strategy = inputControls.strategy;
controls.Nlive = inputControls.Nlive;
controls.nmcmc = inputControls.Nmcmc;
controls.propScale = inputControls.propScale;
controls.nsTolerance = inputControls.nsTolerance;
switch inputControls.calcSldDuringFit
    case 'no'
        controls.calcSld = 0;
    otherwise
        controls.calcSld = 1;
end
controls.repeats = inputControls.repeats;
controls.nsimu = inputControls.nsimu;
controls.burnin = inputControls.burnin;

%Also need to deal with the checks...

checks.params_fitYesNo = inputStruct.paramFitYesNo;
checks.backs_fitYesNo = inputStruct.backParFitYesNo;
checks.shifts_fitYesNo = inputStruct.qzshiftFitYesNo;
checks.scales_fitYesNo = inputStruct.scalefactorFitYesNo;
checks.nbairs_fitYesNo = inputStruct.nbaFitYesNo;
checks.nbsubs_fitYesNo = inputStruct.nbsFitYesNo;
checks.resol_fitYesNo = inputStruct.resolFitYesNo;

controls.checks = checks;

%********* Old API problem structure******
%              resample: [7×1 double]
%     numberOfContrasts: 7
%              geometry: 'Air / Liquid (or solid)'
%                nbairs: [0 0 0 0 0 0 0]
%                nbsubs: [0 6.3500e-06 0 6.3500e-06 0 6.3500e-06 6.3500e-06]
%         contrastBacks: [2 1 2 1 2 1 1]
%        contrastShifts: [1 1 1 1 1 1 1]
%        contrastScales: [1 1 1 1 1 1 1]
%          contrastNbas: [1 1 1 1 1 1 1]
%          contrastNbss: [2 1 2 1 2 1 1]
%           contrastRes: [1 1 1 1 1 1 1]
%                 backs: [2.9718e-06 5.1473e-06]
%                shifts: 0
%                    sf: 0.2325
%                   nba: 0
%                   nbs: [6.3500e-06 0]
%                   res: 0.0300
%           dataPresent: [1 1 1 1 1 1 1]
%               nParams: 10
%                params: [6.9950 18.6911 6.9473e-06 3.0000 -2.2818e-07 7.0421 5.8839e-06 3.0000 1.8126e-06 10.0536]
%        numberOfLayers: 4
%             whichType: 1
%            fileHandle: []
%                  lang: []
%                module: []
%               fitpars: [18×1 double]
%             otherpars: [0×1 double]
%             fitconstr: [18×2 double]
%           otherconstr: [0×2 double]
%             fitConstr: [18×2 double]


%old API limits structure....
% problemDef_limits = 
% 
%   struct with fields:
% 
%     params: [10×2 double]
%      backs: [2×2 double]
%     scales: [0.1000 0.4000]
%     shifts: [-0.0300 0.0300]
%        nba: [-1.0000e-05 1.0000e-05]
%        nbs: [2×2 double]
%        res: [0.0100 0.0500]










%new API style problem struct....









end