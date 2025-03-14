function [qzshifts,scalefactors,bulkIns,bulkOuts,chis,reflectivity,...
    simulation,shiftedData,backgrounds,resolutions,layerSlds,sldProfiles,...
    resampledLayers,subRoughs] = standardLayers(problemStruct,controls)
    % This is the main reflectivity calculation of the standard layers
    % calculation type. It extracts the required parameters for the contrasts
    % from the input arrays, then passes the main calculation to
    % 'coreLayersCalculation', which carries out the calculation itself. 
    % The core calculation is common for both standard and custom layers.

    % Extract parameters from problemStruct
    [numberOfContrasts, geometry, contrastBackgroundIndices, contrastQzshiftIndices,...
     contrastScalefactorIndices, contrastBulkInIndices, contrastBulkOutIndices,...
     contrastResolutionIndices, ~, backgroundParamValues, qzshiftValues,...
     scalefactorValues, bulkInValues, bulkOutValues, resolutionParamValues,...
     ~, dataPresent, nParams, paramValues, ~, resample, contrastBackgroundTypes,...
     contrastBackgroundActions, contrastResolutionTypes, ~, useImaginary,...
     repeatLayers, data, dataLimits, simulationLimits, contrastLayersIndices,...
     layersDetails, customFiles, ~] = extractProblemParams(problemStruct);
    
    calcSld = controls.calcSldDuringFit;
    parallel = controls.parallel;
    resampleMinAngle = controls.resampleMinAngle;
    resampleNPoints = controls.resampleNPoints;
    
    % Allocate the memory for the output arrays before the main loop
    qzshifts = zeros(numberOfContrasts,1);
    scalefactors = zeros(numberOfContrasts,1);
    bulkIns = zeros(numberOfContrasts,1);
    bulkOuts = zeros(numberOfContrasts,1);
    chis = zeros(numberOfContrasts,1);
    subRoughs = zeros(numberOfContrasts,1);
    
    reflectivity = cell(numberOfContrasts,1);  
    simulation = cell(numberOfContrasts,1);
    shiftedData = cell(numberOfContrasts,1);
    backgrounds = cell(numberOfContrasts,1);
    resolutions = cell(numberOfContrasts,1);
    layerSlds = cell(numberOfContrasts,1);
    sldProfiles = cell(numberOfContrasts,1);
    resampledLayers = cell(numberOfContrasts,1);
    
    % First we need to allocate the absolute values of the input
    % parameters to all the layers in the layers list. This only needs
    % to be done once, and so is done outside the contrasts loop
    layerValues = allocateParamsToLayers(paramValues, layersDetails);   
    
    % Substrate roughness is always first parameter for standard layers
    for i = 1:numberOfContrasts
        subRoughs(i) = paramValues(1);
    end
    
    if strcmpi(parallel, coderEnums.parallelOptions.Contrasts)
    
        % Loop over all the contrasts
        parfor i = 1:numberOfContrasts
    
            [qzshifts(i),scalefactors(i),bulkIns(i),bulkOuts(i),chis(i),...
             reflectivity{i},simulation{i},shiftedData{i},backgrounds{i},...
             resolutions{i},layerSlds{i},sldProfiles{i},resampledLayers{i}...
             ] = contrastCalculation(contrastBackgroundIndices{i},...
             contrastQzshiftIndices(i),contrastScalefactorIndices(i),...
             contrastBulkInIndices(i),contrastBulkOutIndices(i),...
             contrastResolutionIndices{i},backgroundParamValues,qzshiftValues,...
             scalefactorValues,bulkInValues,bulkOutValues,resolutionParamValues,...
             dataPresent(i),data{i},dataLimits{i},simulationLimits{i},repeatLayers{i},...
             contrastBackgroundTypes{i},contrastBackgroundActions{i},...
             contrastResolutionTypes{i},customFiles,nParams,parallel,...
             resampleMinAngle,resampleNPoints,resample(i),geometry,...
             subRoughs(i),calcSld,contrastLayersIndices{i},layerValues);
    
        end
        
    else
    
        % Loop over all the contrasts
        for i = 1:numberOfContrasts
            
            [qzshifts(i),scalefactors(i),bulkIns(i),bulkOuts(i),chis(i),...
             reflectivity{i},simulation{i},shiftedData{i},backgrounds{i},...
             resolutions{i},layerSlds{i},sldProfiles{i},resampledLayers{i}...
             ] = contrastCalculation(contrastBackgroundIndices{i},...
             contrastQzshiftIndices(i),contrastScalefactorIndices(i),...
             contrastBulkInIndices(i),contrastBulkOutIndices(i),...
             contrastResolutionIndices{i},backgroundParamValues,qzshiftValues, ...
             scalefactorValues,bulkInValues,bulkOutValues,resolutionParamValues,...
             dataPresent(i),data{i},dataLimits{i},simulationLimits{i},repeatLayers{i},...
             contrastBackgroundTypes{i},contrastBackgroundActions{i},...
             contrastResolutionTypes{i},customFiles,nParams,parallel,...
             resampleMinAngle,resampleNPoints,resample(i),geometry,...
             subRoughs(i),calcSld,contrastLayersIndices{i},layerValues);

        end
    end

    % Remove dummy imaginary column if present
    if ~useImaginary
        for i=1:numberOfContrasts
            layerSlds{i}(:,3) = [];
            resampledLayers{i}(:,3) = [];
        end
    end

end


function [qzshift,scalefactor,bulkIn,bulkOut,chi,reflectivity,simulation,...
    shiftedData,background,resolution,layerSld,sldProfile,...
    resampledLayers] = contrastCalculation(backgroundParamIndex,qzshiftIndex,...
    scalefactorIndex,bulkInIndex,bulkOutIndex,resolutionParamIndex,...
    backgroundParamValues,qzshiftValues,scalefactorValues,bulkInValues,...
    bulkOutValues,resolutionParamValues,dataPresent,data,dataLimits,...
    simulationLimits,repeatLayers,backgroundType,backgroundAction,...
    resolutionType,customFiles,nParams,parallel,resampleMinAngle,resampleNPoints,...
    resample,geometry,roughness,calcSld,layerIndices,layerValues)

    % Extract the relevant parameter values for this contrast
    % from the input arrays.
    % First need to decide which values of the backgrounds, scalefactors
    % data shifts and bulk contrasts are associated with this contrast
    [qzshift,scalefactor,bulkIn,bulkOut] = backSort( ...
        qzshiftIndex,scalefactorIndex,bulkInIndex,bulkOutIndex, ...
        qzshiftValues,scalefactorValues,bulkInValues,bulkOutValues);
    
    % Also need to determine which layers from the overall layers list
    % are required for this contrast, and put them in the correct order 
    % according to geometry
    layers = allocateLayersForContrast(layerIndices,layerValues);

    % Apply scale factors and q shifts to the data
    shiftedData = shiftData(scalefactor,qzshift,dataPresent,data,dataLimits,simulationLimits);
    [simulationXData, dataIndices] = makeSimulationRange(shiftedData, simulationLimits);

    background = constructBackground(backgroundType,backgroundParamIndex, ...
        shiftedData,customFiles,backgroundParamValues,simulationXData,dataIndices);
    resolution = constructResolution(resolutionType,resolutionParamIndex, ...
        shiftedData,customFiles,resolutionParamValues,simulationXData,dataIndices);

    % Call the core layers calculation
    [reflectivity,simulation,shiftedData,layerSld,sldProfile,...
        resampledLayers] = normalTF.coreLayersCalculation(layers,roughness,...
        geometry,bulkIn,bulkOut,resample,calcSld,shiftedData,simulationXData,dataIndices,repeatLayers,...
        resolution,background,backgroundAction,parallel,resampleMinAngle,resampleNPoints);

    % Calculate chi squared
    chi = chiSquared(shiftedData,reflectivity,nParams);

end
