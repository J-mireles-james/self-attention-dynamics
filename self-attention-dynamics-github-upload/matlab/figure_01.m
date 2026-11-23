% Figure 1, Section 3.1 — periodic orbit and selected trajectories.
% Author: J.D. Mireles James
% AI use: Created and tested with AI assistance. See README for test scope.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'functions'));
dataDir = fullfile(root,'data');
recompute = false; % Set true to repeat the numerical experiment.
load(fullfile(dataDir,'ex5_periodic_portrait_sync.mat'),'data');
portrait(data,recompute);
