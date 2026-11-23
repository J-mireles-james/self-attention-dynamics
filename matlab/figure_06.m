% Figure 6, Section 3.3. Forward integration scan.
% Author: J.D. Mireles James
% AI use: Created and tested with AI assistance. See README for test scope.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'functions'));
dataDir = fullfile(root,'data');
recompute = false; % Set true to repeat the numerical experiment.
figure;
parameter_scan(fullfile(dataDir,'ex5_complex_periodic_bifurcation_left.mat'),'leftResult',recompute);
parameter_scan(fullfile(dataDir,'ex5_complex_periodic_bifurcation_right.mat'),'rightResult',recompute);
