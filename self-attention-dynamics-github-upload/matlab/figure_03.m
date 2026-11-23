% Figure 3, Section 3.2. Fixed-parameter diagnostics at beta = 0.1.
% Author: J.D. Mireles James
% AI use: Created and tested with AI assistance. See README for test scope.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'functions'));
dataDir = fullfile(root,'data');
recompute = false; % Set true to repeat the numerical experiment.
diagnostics(fullfile(dataDir,'ex5_fixed_beta_0p1_data.mat'),recompute);
