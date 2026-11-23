% Table 2. Full-duration integration into a fresh temporary folder.
% Author: J.D. Mireles James
% AI use: Created and tested with AI assistance. See README for test scope.

function run_60_particle_heads_definitive_safe(resultFolder)
%% USER SETTINGS
runMode='computeOnly';          % 'computeOnly' or 'reportOnly'
if ~exist(resultFolder,'dir'),mkdir(resultFolder);end

% Both original independent seeds are available. Use [0] for the first
% definitive run, then [1] for independent confirmation.
seedVariants = [0 1];
includeIsolatedBaseline = true;
includeMasked60 = false; % Optional additional full-length masked control

autoResume=true;
forceRecompute=false;

% Refined h and extended times, significantly more demanding than
% the original Python h=.025, burn=1000, Lyapunov-time=4000 test.
cfg.h=0.0125;              % half the original step (RK4)
cfg.transientTime=2000;    % original 1000
cfg.observationTime=16000; % original 4000
cfg.qrTime=1;
cfg.historyTime=40;
cfg.checkpointTime=250;
cfg.runMode=runMode;
cfg.resumeFromCheckpoint=autoResume;
cfg.forceRecompute=forceRecompute;
cfg.checkJacobian=true;
cfg.outputFolder=resultFolder;
cfg.summaryPrefix='sixty_particle_table';

Astar=[-0.05,-2.30;2.30,-0.05];
Bstar=[-1.20,0.40;-0.80,-2.40];
A=repmat(Astar,1,1,20);
B=repmat(Bstar,1,1,20);
beta=1.25*ones(1,20);
eps=0.1;

jobs = struct('name',{},'topology',{},'epsilon',{}, ...
              'A',{},'B',{},'beta',{},'theta0',{},'k',{});
for variant=seedVariants
    theta0=attentionLE_reference_seeds(sprintf('sixty_v%d',variant));
    job=struct('name',sprintf('M20_heads_eps0p1_seed%d',variant), ...
        'topology','heads','epsilon',eps,'A',A,'B',B, ...
        'beta',beta,'theta0',theta0,'k',40);
    jobs(end+1)=job; %#ok<SAGROW>
    if includeMasked60
        job.topology='mask';
        job.name=sprintf('M20_mask_eps0p1_seed%d',variant);
        jobs(end+1)=job; %#ok<SAGROW>
    end
end

if includeIsolatedBaseline
    % The uncoupled M=20 product has 20 copies of the isolated 3-particle
    % spectrum. This run measures ONE such subsystem, not all 60 LEs.
    theta0=attentionLE_reference_seeds('sixty_v0');
    job=struct('name','M1_uncoupled_single_group_baseline', ...
        'topology','mask','epsilon',0,'A',Astar,'B',Bstar, ...
        'beta',1.25,'theta0',theta0(1:3),'k',3);
    jobs(end+1)=job;
end

attentionLE_core(cfg,jobs);

% The original exploratory 60-particle result at h=.025, burn=1000,
% obs=4000 had 26 exponents > 0.01 from each of two distinct seeds.
% Refined values may differ, especially close to zero or across basins.
