% Figure 2, Section 3.2. Forward integration scan.
% Author: J.D. Mireles James
% AI use: Created and tested with AI assistance. See README for test scope.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'functions'));
dataDir = fullfile(root,'data');
recompute = false; % Set true to repeat the numerical experiment.
figure;
parameter_scan(fullfile(dataDir,'ex5_beta0_full_data.mat'),'fullResult',recompute);
