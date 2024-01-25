function out = orsoTest1()
% ORSO validation Test1

layers = dlmread('test1.layers');

% Change the units to Å
layers(:,2) = layers(:,2) .* 1e-6;
layers(:,3) = layers(:,3) .* 1e-6;

% Read in the data.....
data = dlmread('test1.dat');

% Group the Layers
thick = layers(:,1);
sld = complex(layers(:,2),layers(:,3));
rough = layers(:,4);

% Calculate reflectivity...
q = data(:,1);
N = size(layers,1);
ref = abelesSingle(q,N,thick,sld,rough);

% Plot the comparison....
figure(1); clf
semilogy(q,ref,'k-','LineWidth',2)
hold on
plot(data(:,1),data(:,2),'r.')

% Calculate the output....
out = sum(sum((data(:,2) - ref).^2));

