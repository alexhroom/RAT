% Test a simple example using simulated data.
% The 'data' in this example is actually simulation data using true values:
% Layer |  Thickness | SLD  | Roughness 
% ------|------------|------|-----------
%   1   |     60     | 4e-6 |     6     
%   2   |     150    | 1e-6 |     8     
%

function results = runtest()
problem = createProject(name='optim test');

problem.showPriors = true;

% true values in order: thickness, SLD, roughness for layer 1 then 2
param_values = [60 4e-6 2; 150 1e-6 8];

perturb = true;
max_ptb = 0.25;
if perturb
    % perturb parameter values by some random percentage from 0% to max_ptb%
    ptb = rand(1, 6) * max_ptb;
    signs = randi([0 1], 1, 6);
    for i = 1:6
        if signs(i)
            param_values(i) = param_values(i) + param_values(i)*ptb(i);
        else
            param_values(i) = param_values(i) - param_values(i)*ptb(i);
        end
    end
end

paramGroup = {
               {'Layer 1 Thickness', 40, param_values(1, 1), 80, true, 'gaussian', 60, 60*max_ptb};
               {'Layer 1 SLD', 3e-6, param_values(1, 2), 5e-6, true, 'uniform', 0, Inf};
               {'Layer 1 Roughness', 1, param_values(1, 3), 5, true, 'uniform', 0, Inf};
               {'Layer 2 Thickness', 100, param_values(2, 1), 200, true, 'gaussian', 150, 150*max_ptb};
               {'Layer 2 SLD', 5e-7, param_values(2, 2), 1.5e-6, true, 'uniform', 0, Inf};
               {'Layer 2 Roughness', 6, param_values(2, 3), 10, true, 'uniform', 0, Inf};
             };

problem.addParameterGroup(paramGroup);

layers = {
    {'Layer 1', 'Layer 1 Thickness', 'Layer 1 SLD', 'Layer 1 Roughness'};
    {'Layer 2', 'Layer 2 Thickness', 'Layer 2 SLD', 'Layer 2 Roughness'};
};

problem.addLayerGroup(layers);

stack = {'Layer 1', 'Layer 2'};

data_array = readmatrix("data.dat");
problem.addData("dataset", data_array);

problem.addContrast( ...
    'name',         'Contrast 1', ...
    'BulkIn',       'SLD Air',...
    'BulkOut',      'SLD D2O',...
    'background',   'Background 1',...
    'resolution',   'Resolution 1',... 
    'scalefactor',  'Scalefactor 1',...
    'data',         'dataset',...
    'model',        stack);


problem.setScalefactor(1, "max", 1);
problem.setScalefactor(1, "value", 1);

procedures = ["calculate", "simplex", "de", "ns", "dream"];
results = zeros(1, 5);
parfor i=1:5
controls = controlsClass();
controls.procedure = procedures(i);
[~, r] = RAT(problem, controls);
results(i) = r.calculationResults.sumChi;
end
end

chis = zeros(10, 5);
for repeat=1:10
    chis(repeat, :) = runtest();
end