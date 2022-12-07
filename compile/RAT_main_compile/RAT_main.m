function [outProblemDef,problem,results,bayesResults] = RAT_main(problemDef,problemDef_cells,problemDef_limits,controls,priors)


result = cell(6,1);
for i = 1:6
    result{i} = {[0]};
end
results = result;

numberOfContrasts = problemDef.numberOfContrasts;
preAlloc = zeros(numberOfContrasts,1);

problem = struct('ssubs',preAlloc,...
    'backgrounds',preAlloc,...
    'qshifts',preAlloc,...
    'scalefactors',preAlloc,...
    'nbairs',preAlloc,...
    'nbsubs',preAlloc,...
    'resolutions',preAlloc,...
    'calculations',struct('all_chis',preAlloc,'sum_chi',0),...
    'allSubRough',preAlloc,...
    'resample',preAlloc);

% Make empty bayes results even though we may not fill it (for output purposes)
bayesResults.res = [];
bayesResults.chain = [];
bayesResults.s2chain = [];
bayesResults.ssChain = [];
bayesResults.bestPars = [];
outProblemDef = problemDef;

%Decide what we are doing....
action = controls.proc;
switch lower(action)
    case 'calculate' %Just a single reflectivity calculation
        [problem,results] = singleCalculation(problemDef,problemDef_cells,problemDef_limits,controls);
        outProblemDef = problemDef;
    case 'simplex'
        if ~strcmpi(controls.display,'off')
            ratSendTextOutput(sprintf('\nRunning simplex\n\n'));
        end
        [outProblemDef,problem,results] = runSimplex(problemDef,problemDef_cells,problemDef_limits,controls);
    case 'de'
        if ~strcmpi(controls.display,'off')
            ratSendTextOutput(sprintf('\nRunning Differential Evolution\n\n'));
        end
        [outProblemDef,problem,results] = runDE(problemDef,problemDef_cells,problemDef_limits,controls);
        %     case 'bayes'
        %         if ~strcmpi(controls.display,'off')
        %             ratSendTextOutput(sprintf('\nRunning DRAM\n\n'));
        %         end
        %         [outProblemDef,problem,results,bayesResults] = runDram(problemDef,problemDef_cells,problemDef_limits,controls,priors);
        %     case 'ns'
        %         if ~strcmpi(controls.display,'off')
        %             ratSendTextOutput(sprintf('\nRunning Nested Sampler\n\n'));
        %         end
        %         [outProblemDef,problem,results,bayesResults] = runNestedSampler(problemDef,problemDef_cells,problemDef_limits,controls);
        %     case 'paramonte'
        %         % TODO
end

end