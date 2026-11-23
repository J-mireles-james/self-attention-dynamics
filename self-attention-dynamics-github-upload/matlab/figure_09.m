% Figure 9, Section 5.1. Six-particle Lyapunov spectra.
% Author: J.D. Mireles James
% AI use: Created and tested with AI assistance. See README for test scope.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'functions'));
dataDir = fullfile(root,'data');
recompute = false; % Set true to repeat the numerical experiment.
figure; tiledlayout(3,2);
kinds={'mask','heads'}; tags={'0p001','0p1'};
for row=1:3
    for col=1:2
        if row~=2
            kind=kinds{1+(row==3)};
            file=sprintf('ex5_coupled_six_long_%s_eps_%s.mat',kind,tags{col});
            S=load(fullfile(dataDir,file)); R=S.R; c=S.configSaved;
            if recompute
                R=coupled_compute_case(c,R.epsilon,kind,round(c.transientTime/c.h), ...
                    round(c.observationTime/c.qrTime),round(c.qrTime/c.h), ...
                    round(c.historyTime/c.qrTime),250,'',false);
            end
            if row==1, masked{col}=R; end
        else
            % Equation (15) identifies the label bias with the block mask.
            R=masked{col};
        end
        nexttile; plot(R.times,R.history.'); xlabel('Time'); ylabel('Finite-time exponents');
        names={'Block mask','Group labels','Separate heads'};
        title(sprintf('%s, epsilon = %g',names{row},R.epsilon));
    end
end
