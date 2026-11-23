% Sections 3.2–3.3. Warm-start scans using signed section crossings.
% Author: J.D. Mireles James
% AI use: Created and tested with AI assistance. See README for test scope.

function parameter_scan(file,key,recompute)
S = load(file,key,'params');
r = S.(key);
p = S.params;
if recompute
    if isfield(p,'seedTheta'), seed=p.seedTheta; else, seed=p.initialTheta; end
    opts=odeset('RelTol',p.relativeTolerance,'AbsTol',p.absoluteTolerance,'MaxStep',p.maxStep);
    events=odeset(opts,'Events',@local_section_event);
    settings.burn=r.transientTimeUsed;
    settings.observationTime=r.observationTime;
    r=local_scan(r.betas,seed,p.A,p.B,@ode113,opts,events,settings);
end
points=r.downPoints;
plot(points(:,1),points(:,3),'.','MarkerSize',1,'Color',[0 .25 .55]);
hold on; xlabel('\beta'); ylabel('\theta_3'); ylim([0 2*pi]);
end
function result=local_scan(betas,x0,A,B,solver,opts,optsEvents,settings)
% Warm-start forward flow in the specified order, from the central seed.
    n=numel(betas);
    downCells=cell(n,1);upCells=cell(n,1);
    terminalStates=nan(3,n);
    speeds=nan(n,1);
    nDown=zeros(n,1);nUp=zeros(n,1);
    burnUsed=zeros(n,1);
    x=mod(x0(:),2*pi);
    reportEvery=max(1,round(n/15));

    for j=1:n
        beta=betas(j);
        burn=settings.burn(j);
        burnUsed(j)=burn;
        rhs=@(~,theta)local_field(theta,A,B,beta);
        % No event detection in the discarded transient.
        [~,yBurn]=solver(rhs,[0,burn],x,opts);
        xBurn=mod(yBurn(end,:).',2*pi);

        % MATLAB Events locate theta1=0 (mod 2*pi); classify physical
        % crossing orientation via sign of F1 rather than event direction.
        [~,yObs,tEvents,yEvents]=solver(rhs,[0,settings.observationTime], ...
                                     xBurn,optsEvents);
        x=mod(yObs(end,:).',2*pi);
        terminalStates(:,j)=x;
        speeds(j)=norm(local_field(x,A,B,beta),inf);
        if isempty(tEvents)
            down=zeros(0,2);up=zeros(0,2);
        else
            good=tEvents>1e-8;
            yEvents=yEvents(good,:);
            crossings=zeros(size(yEvents,1),1);
            for k=1:size(yEvents,1)
                v=local_field(yEvents(k,:).',A,B,beta);
                crossings(k)=v(1);
            end
            down=mod(yEvents(crossings<-1e-9,2:3),2*pi);
            up=mod(yEvents(crossings>1e-9,2:3),2*pi);
        end
        nDown(j)=size(down,1);nUp(j)=size(up,1);
        downCells{j}=[repmat(beta,nDown(j),1),down];
        upCells{j}=[repmat(beta,nUp(j),1),up];
        if j==1 || j==n || mod(j,reportEvery)==0
            fprintf('  %4d/%4d beta=%.6f down=%3d up=%3d ||F(end)||=%.2e\n', ...
                j,n,beta,nDown(j),nUp(j),speeds(j));
        end
    end
    result=struct('betas',betas(:), ...
        'downPoints',vertcat(downCells{:}), ...
        'upPoints',vertcat(upCells{:}), ...
        'nDown',nDown,'nUp',nUp, ...
        'terminalStates',terminalStates, ...
        'terminalFieldNormInf',speeds, ...
        'transientTimeUsed',burnUsed, ...
        'observationTime',settings.observationTime, ...
        'startTheta',mod(x0(:),2*pi));
end

function f=local_field(theta,A,B,beta)
% Original Ex5 model; vectorized across three particles.
    theta=theta(:);
    u=[cos(theta).';sin(theta).'];
    v=[-u(2,:);u(1,:)];
    scores=beta*(u.'*B*u);
    scores=scores-max(scores,[],2);
    w=exp(scores);
    w=w./sum(w,2);
    meanVectors=u*w.';
    f=sum(v.*(A*meanVectors),1).';
end

function [value,isterminal,direction]=local_section_event(~,theta)
% theta1 = 2*pi*k; real orientation is computed afterward from F1.
    value=sin(theta(1)/2);
    isterminal=0;
    direction=0;
end

