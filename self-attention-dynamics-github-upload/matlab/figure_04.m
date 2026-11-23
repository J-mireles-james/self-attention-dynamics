% Figure 4, Section 3.2. Fixed-parameter diagnostics at beta = 1.25.
% Author: J.D. Mireles James
% AI use: Created and tested with AI assistance. See README for test scope.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'functions'));
dataDir = fullfile(root,'data');
recompute = false; % Set true to repeat the numerical experiment.
diagnostics(fullfile(dataDir,'ex5_fixed_beta_1p25_data.mat'),recompute);
