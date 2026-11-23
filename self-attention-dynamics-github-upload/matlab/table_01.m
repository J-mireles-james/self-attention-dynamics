% Table 1, Section 5.2. Leading exponents and counts above 0.01.
% Author: J.D. Mireles James
% AI use: Created and tested with AI assistance. See README for test scope.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'functions'));
dataDir = fullfile(root,'data');
recompute = false; % Set true to repeat the numerical experiment.
if recompute
    folder=tempname; mkdir(folder);
    run_9_particle_mixed_safe(folder);
else
    folder=dataDir;
end
T=readtable(fullfile(folder,'nine_particle_table.csv'),'VariableNamingRule','preserve');
disp(T(:,{'case','particles','epsilon','k','lambda_01','lambda_02','positive_gt_0p01'}));
