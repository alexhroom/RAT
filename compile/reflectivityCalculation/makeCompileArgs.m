function ARGS = makeCompileArgs()

% Define the arguments for compiling reflectivityCalculation
% using codegen.

%% Define argument types for entry-point 'reflectivityCalculation'.
maxArraySize = 10000;
maxDataSize = 10000;

ARGS = cell(1,1);
ARGS{1} = cell(3,1);
ARGS_1_1 = struct;
ARGS_1_1.TF = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_1.resample = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.dataPresent = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.oilChiDataPresent = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.numberOfContrasts = coder.typeof(0);
ARGS_1_1.geometry = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_1.useImaginary = coder.typeof(true,[1 1],[0 0]);
ARGS_1_1.contrastBackgroundParams = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.contrastBackgroundActions = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.contrastQzshifts = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.contrastScalefactors = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.contrastBulkIns = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.contrastBulkOuts = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.contrastResolutionParams = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.backgroundParams = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.qzshifts = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.scalefactors = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.bulkIn = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.bulkOut = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.resolutionParams = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.params = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.numberOfLayers = coder.typeof(0);
ARGS_1_1.modelType = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_1.contrastCustomFiles = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.contrastDomainRatios = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.domainRatio = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_1.numberOfDomainContrasts = coder.typeof(0);
ARGS_1_1.fitParams = coder.typeof(0,[maxArraySize 1],[1 0]);
ARGS_1_1.otherParams = coder.typeof(0,[maxArraySize 1],[1 0]);
ARGS_1_1.fitLimits = coder.typeof(0,[maxArraySize 2],[1 0]);
ARGS_1_1.otherLimits = coder.typeof(0,[maxArraySize 2],[1 0]);
ARGS{1}{1} = coder.typeof(ARGS_1_1);
ARGS_1_2 = cell([1 21]);
ARG = coder.typeof(0,[1 2]);
ARGS_1_2{1} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof(0,[maxDataSize  5],[1 1]);
ARGS_1_2{2} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof(0,[1 2]);
ARGS_1_2{3} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof(0,[1 2]);
ARGS_1_2{4} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof(0,[1 maxArraySize],[1 1]);
ARGS_1_2{5} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof(0,[1 10],[1 1]);
ARGS_1_2{6} = coder.typeof({ARG}, [maxArraySize 1],[1 0]);
ARG = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_2{7} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_2{8} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_2{9} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_2{10} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_2{11} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_2{12} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_2{13} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_2{14} = coder.typeof({ARG}, [1 maxArraySize], [0 1]);
ARG = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_2{15} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_2{16} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof(0,[maxArraySize  5],[1 1]);
ARGS_1_2{17} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof(0,[1 2]);
ARGS_1_2{18} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof(0,[1 maxArraySize],[1 1]);
ARGS_1_2{19} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_2{20} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARG = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_2{21} = coder.typeof({ARG}, [1 maxArraySize],[0 1]);
ARGS{1}{2} = coder.typeof(ARGS_1_2,[1 21]);
ARGS{1}{2} = ARGS{1}{2}.makeHeterogeneous();
ARGS_1_3 = struct;
ARGS_1_3.param = coder.typeof(0,[maxArraySize 2],[1 0]);
ARGS_1_3.backgroundParam = coder.typeof(0,[maxArraySize 2],[1 0]);
ARGS_1_3.scalefactor = coder.typeof(0,[maxArraySize 2],[1 0]);
ARGS_1_3.qzshift = coder.typeof(0,[maxArraySize 2],[1 0]);
ARGS_1_3.bulkIn = coder.typeof(0,[maxArraySize 2],[1 0]);
ARGS_1_3.bulkOut = coder.typeof(0,[maxArraySize 2],[1 0]);
ARGS_1_3.resolutionParam = coder.typeof(0,[maxArraySize 2],[1 0]);
ARGS_1_3.domainRatio = coder.typeof(0,[maxArraySize 2],[1 0]);
ARGS{1}{3} = coder.typeof(ARGS_1_3);
ARGS_1_4 = struct;
ARGS_1_4.procedure = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_4.parallel = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_4.resampleMinAngle = coder.typeof(0);
ARGS_1_4.resampleNPoints = coder.typeof(0);
ARGS_1_4.calcSldDuringFit = coder.typeof(true);
ARGS_1_4.display = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_4.xTolerance = coder.typeof(0);
ARGS_1_4.funcTolerance = coder.typeof(0);
ARGS_1_4.maxFuncEvals = coder.typeof(0);
ARGS_1_4.maxIterations = coder.typeof(0);
ARGS_1_4.updateFreq = coder.typeof(0,[1 1]);
ARGS_1_4.updatePlotFreq = coder.typeof(0,[1 1]);
ARGS_1_4.populationSize = coder.typeof(0);
ARGS_1_4.fWeight = coder.typeof(0);
ARGS_1_4.crossoverProbability = coder.typeof(0);
ARGS_1_4.strategy = coder.typeof(0);
ARGS_1_4.targetValue = coder.typeof(0);
ARGS_1_4.numGenerations = coder.typeof(0);
ARGS_1_4.nLive = coder.typeof(0);
ARGS_1_4.nMCMC = coder.typeof(0);
ARGS_1_4.propScale = coder.typeof(0);
ARGS_1_4.nsTolerance = coder.typeof(0);
ARGS_1_4.nSamples = coder.typeof(0,[1 1]);
ARGS_1_4.nChains = coder.typeof(0,[1 1]);
ARGS_1_4.jumpProbability = coder.typeof(0,[1 1]);
ARGS_1_4.pUnitGamma = coder.typeof(0,[1 1]);
ARGS_1_4.boundHandling = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS_1_4.adaptPCR = coder.typeof(true);
ARGS_1_4_checks = struct;
ARGS_1_4_checks.fitParam = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_4_checks.fitBackgroundParam = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_4_checks.fitQzshift = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_4_checks.fitScalefactor = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_4_checks.fitBulkIn = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_4_checks.fitBulkOut = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_4_checks.fitResolutionParam = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_4_checks.fitDomainRatio = coder.typeof(0,[1 maxArraySize],[0 1]);
ARGS_1_4.checks = coder.typeof(ARGS_1_4_checks);
ARGS_1_4.IPCFilePath = coder.typeof('X',[1 maxArraySize],[0 1]);
ARGS{1}{4} = coder.typeof(ARGS_1_4);

end
