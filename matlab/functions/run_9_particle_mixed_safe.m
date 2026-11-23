% Table 1. Full-duration integration into a fresh temporary folder.
% Author: J.D. Mireles James
% AI use: Created and tested with AI assistance. See README for test scope.

function run_9_particle_mixed_safe(resultFolder)
%% USER SETTINGS
runMode = 'computeOnly'; % 'computeOnly' or 'reportOnly'
if ~exist(resultFolder,'dir'), mkdir(resultFolder); end

% Use [0] to reproduce just the original first initial state, or [0 1]
% for the two independent initial configurations used in the paper study.
seedVariants = [0 1];
epsilonList = 0.1; % Original nine-particle experiment. Optional [0.001 0.1].

autoResume = true;
forceRecompute = false; % true removes previous results/checkpoints for jobs

% Original verified long-run parameters:
cfg.h=0.02;
cfg.transientTime=1500;
cfg.observationTime=9000;
cfg.qrTime=1;
cfg.historyTime=20;
cfg.checkpointTime=250;
cfg.runMode=runMode;
cfg.resumeFromCheckpoint=autoResume;
cfg.forceRecompute=forceRecompute;
cfg.checkJacobian=true;
cfg.outputFolder=resultFolder;
cfg.summaryPrefix='nine_particle_table';

Astar = [-0.05, -2.30; 2.30, -0.05];
Bstar = [-1.20, 0.40; -0.80, -2.40];
A = repmat(Astar,1,1,3);
B = repmat(Bstar,1,1,3);
betas = [1.25, 0.55, 0.55];

jobs = struct('name',{},'topology',{},'epsilon',{}, ...
              'A',{},'B',{},'beta',{},'theta0',{},'k',{});
for variant=seedVariants
    theta0=attentionLE_reference_seeds(sprintf('nine_v%d',variant));
    for eps=epsilonList
        for topology={'mask','heads'}
            topo=topology{1};
            epTag=strrep(sprintf('%.8g',eps),'.','p');
            job=struct('name',sprintf('C_2P_%s_eps%s_seed%d',topo,epTag,variant), ...
                'topology',topo,'epsilon',eps,'A',A,'B',B, ...
                'beta',betas,'theta0',theta0,'k',9);
            jobs(end+1)=job; %#ok<SAGROW>
        end
    end
end

attentionLE_core(cfg,jobs);

% Expected *qualitative* result at eps=.1 from original Python (NOT
% hard-coded output): mask -> no clearly positive LE; heads -> 3 clearly
% positive LEs (~.14, .095, .042). Long transients and basins matter.
