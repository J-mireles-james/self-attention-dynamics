% Figures 3, 4, 7 and 8. Variational flow, oriented section and angular orbit.
% Author: J.D. Mireles James
% AI use: Created and tested with AI assistance. See README for test scope.

function diagnostics(file,recompute)
S = load(file,'params','result');
p = S.params;
r = S.result;
if recompute
    r = fixed_simulate(p.initialTheta,p.A,p.B,p.beta,@ode113,p,[]);
end
figure; tiledlayout(1,3);
for kind = {'lyapunov','poincare','phase'}
    fixed_plot_one(r,p,'both',60000,Inf,kind{1},0.55,[0.12 0.32 0.60],35,25);
end
end
function result = fixed_simulate(x,A,B,beta,solver,p,logFile)
    P = 2*pi;
    opts = odeset('RelTol',p.relativeTolerance,'AbsTol',p.absoluteTolerance, ...
        'MaxStep',p.maxStep);
    eventOpts = odeset(opts,'Events',@fixed_section_event);
    rhs = @(~,y) fixed_field(y,A,B,beta);
    rhsAug = @(~,z) fixed_augmented_field(z,A,B,beta);
    fixed_log(logFile,'Removing transient (no events or graphics) ...');
    elapsed = 0;
    lastReport = 0;
    while elapsed<p.transientTime
        dt = min(p.burnChunkTime,p.transientTime-elapsed);
        [tb,Y] = solver(rhs,[0,dt],x,opts);
        assert(abs(tb(end)-dt)<=1e-9*max(1,dt) && all(isfinite(Y(end,:))), ...
            'ODE integration failed during transient at t=%.8g.',elapsed);
        x = mod(Y(end,:).',P);
        elapsed = elapsed+dt;
        if elapsed-lastReport>=max(p.transientTime/10,p.burnChunkTime) || ...
                elapsed>=p.transientTime
            fixed_log(logFile,'Transient %.1f / %.1f',elapsed,p.transientTime);
            lastReport=elapsed;
        end
    end
    postTransientTheta = x;
    fixed_log(logFile,'Transient complete. Integrating orbit and variational equations ...');

    % Each QR block starts from a wrapped angle vector and an orthonormal
    % tangent basis. The ODE is integrated on a continuous lift WITHIN
    % the block; wrapping only between blocks does not affect DF or QR.
    nBlocks = ceil(p.observationTime/p.qrTime);
    edges = (0:nBlocks)*p.qrTime;
    edges(end) = p.observationTime;
    % Protect against a floating-point duplicate last endpoint.
    edges = unique(edges);
    nBlocks = numel(edges)-1;
    Q = eye(3);
    logGrowth = zeros(3,1);
    divergenceIntegral = 0;
    LEhistory = zeros(3,nBlocks);
    LEtimes = edges(2:end);
    meanDivergence = zeros(1,nBlocks);
    cloudCells = cell(nBlocks,1);
    timeCells = cell(nBlocks,1);
    downCells = cell(nBlocks,1);
    upCells = cell(nBlocks,1);
    reportEvery = max(1,ceil(nBlocks/12));
    for k=1:nBlocks
        startTime = edges(k);
        dt = edges(k+1)-startTime;
        % Fixed output sampling, adaptive internal ODE steps. Include
        % at least three output times to request an explicit sample grid.
        tspan = linspace(0,dt,max(2,ceil(dt/p.sampleTime))+1);
        z0 = [x;Q(:);0];  % last component integrates div F on this block
        [ts,Z,te,ZE] = solver(rhsAug,tspan,z0,eventOpts);
        assert(abs(ts(end)-dt)<=1e-9*max(1,dt) && all(isfinite(Z(end,:))), ...
            'ODE integration failed to finish block %d.',k);
        x = mod(Z(end,1:3).',P);
        Phi = reshape(Z(end,4:12),3,3);
        [Q,R] = qr(Phi,0);
        diagonal = abs(diag(R));
        assert(all(diagonal>0 & isfinite(diagonal)), ...
            'Singular QR factor: decrease qrTime or tighten tolerances.');
        logGrowth = logGrowth+log(diagonal);
        divergenceIntegral = divergenceIntegral+Z(end,13);
        LEhistory(:,k) = logGrowth/LEtimes(k);
        meanDivergence(k) = divergenceIntegral/LEtimes(k);
        cloudCells{k} = mod(Z(2:end,1:3).',P);
        timeCells{k} = startTime+ts(2:end).';

        % sin(theta1/2) vanishes at theta1=2*pi*n. Its derivative has
        % alternating sign, so classify using F_1 at the event instead.
        % Exclude the initial event of each block to avoid double counting.
        D = zeros(0,3); U = zeros(0,3);
        if ~isempty(te)
            good = te>1e-8;
            te = te(good);
            XE = ZE(good,1:3);
            speed = zeros(numel(te),1);
            for j=1:numel(te)
                fj = fixed_field(XE(j,:).',A,B,beta);
                speed(j)=fj(1);
            end
            events = [startTime+te(:),mod(XE(:,2:3),P)];
            D = events(speed < -p.crossingSpeedTol,:);
            U = events(speed >  p.crossingSpeedTol,:);
        end
        downCells{k}=D; upCells{k}=U;
        if k==1 || mod(k,reportEvery)==0 || k==nBlocks
            current = sort(LEhistory(:,k),'descend');
            fixed_log(logFile,'t = %8.1f / %.1f   LE = [%+.5f %+.5f %+.5f]', ...
                LEtimes(k),p.observationTime,current);
        end
    end
    LE = sort(LEhistory(:,end),'descend');
    % General Kaplan-Yorke rule, applied to numerical estimates only.
    partial = cumsum(LE);
    j = find(partial>=0,1,'last');
    if isempty(j), KY=0;
    elseif j==3, KY=3;
    elseif LE(j+1)<0, KY=j+partial(j)/abs(LE(j+1));
    else, KY=NaN;
    end
    result = struct();
    result.LE = LE;
    result.LEhistory = LEhistory;
    result.LEtimes = LEtimes;
    result.meanDivergence = meanDivergence;
    result.KYdimension = KY;
    result.sectionDown = vertcat(downCells{:});
    result.sectionUp = vertcat(upCells{:});
    result.cloud = [postTransientTheta,horzcat(cloudCells{:})];
    result.cloudTimes = [0,horzcat(timeCells{:})];
    result.postTransientTheta = postTransientTheta;
    result.terminalTheta = x;
    result.terminalFieldNormInf = norm(fixed_field(x,A,B,beta),inf);
    result.stationaryFlag = isempty(result.sectionDown) && ...
        isempty(result.sectionUp) && result.terminalFieldNormInf<p.stationaryTolerance;
end

function [f,Jac] = fixed_field(theta,A,B,beta)
    % u_i = (cos(theta_i),sin(theta_i)), v_i = u_i'.
    % f_i = v_i' A sum_j alpha_ij u_j; analytic Jacobian as in Ex5.
    theta=theta(:);
    u=[cos(theta).';sin(theta).'];
    v=[-u(2,:);u(1,:)];
    scores=beta*(u.'*B*u);
    weights=exp(scores-max(scores,[],2));
    weights=weights./sum(weights,2);
    means=u*weights.';
    H=A*means;
    f=sum(v.*H,1).';
    if nargout<2, return; end
    Jac=zeros(3);
    for i=1:3
        for ell=1:3
            dm=weights(i,ell)*v(:,ell);
            key=beta*weights(i,ell)*(u(:,i).'*B*v(:,ell));
            dm=dm+key*(u(:,ell)-means(:,i));
            if ell==i
                query=beta*weights(i,:).*(v(:,i).'*B*u);
                dm=dm+(u-means(:,i))*query.';
            end
            Jac(i,ell)=v(:,i).'*A*dm;
            if ell==i
                Jac(i,ell)=Jac(i,ell)-u(:,i).'*H(:,i);
            end
        end
    end
end

function dz = fixed_augmented_field(z,A,B,beta)
    [f,Jac]=fixed_field(z(1:3),A,B,beta);
    Phi=reshape(z(4:12),3,3);
    dPhi=Jac*Phi;
    dz=[f;dPhi(:);trace(Jac)];
end

function err = fixed_check_jacobian(x,A,B,beta)
    [~,Jac]=fixed_field(x,A,B,beta);
    Jfd=zeros(3); h=1e-6;
    for j=1:3
        e=zeros(3,1); e(j)=h;
        Jfd(:,j)=(fixed_field(x+e,A,B,beta)-fixed_field(x-e,A,B,beta))/(2*h);
    end
    err=norm(Jac-Jfd,'fro')/max(1,norm(Jac,'fro'));
end

function [value,isterminal,direction] = fixed_section_event(~,z)
    value=sin(z(1)/2);
    isterminal=0;
    direction=0;
end

function fig = fixed_plot_one(r,p,sectionDirection,maxPhase,maxSection, ...
        kind,phaseWidth,phaseColor,viewAzimuth,viewElevation)
% All graphics are basic MATLAB lines and point markers; no surfaces/lighting.
    P=2*pi;
    ticks=[0,pi/2,pi,3*pi/2,2*pi];
    labels={'0','\pi/2','\pi','3\pi/2','2\pi'};
    tag=sprintf('\\beta = %.8g',p.beta);
    switch lower(kind)
        case 'lyapunov'
            fig=gcf; ax=nexttile;
            h=plot(ax,r.LEtimes,r.LEhistory.','LineWidth',1.1);
            hold(ax,'on'); yline(ax,0,'k:','HandleVisibility','off');
            grid(ax,'on'); box(ax,'on');
            xlabel(ax,'Time after transient');
            ylabel(ax,'Finite-time Lyapunov estimate');
            title(ax,['Lyapunov exponent convergence, ' tag]);
            legend(ax,h,{'QR exponent 1','QR exponent 2','QR exponent 3'}, ...
                'Location','northeast','AutoUpdate','off');
            set(ax,'FontSize',12);
        case 'poincare'
            fig=gcf; ax=nexttile;hold(ax,'on');
            nShown=0;
            if any(strcmpi(sectionDirection,{'down','both'})) && ...
                    ~isempty(r.sectionDown)
                D=r.sectionDown;
                ix=fixed_plot_indices(size(D,1),maxSection);
                plot(ax,D(ix,2),D(ix,3),'.', ...
                    'Color',[0.12,0.32,0.60],'MarkerSize',5, ...
                    'DisplayName','Downward');
                nShown=nShown+numel(ix);
            end
            if any(strcmpi(sectionDirection,{'up','both'})) && ...
                    ~isempty(r.sectionUp)
                U=r.sectionUp;
                ix=fixed_plot_indices(size(U,1),maxSection);
                plot(ax,U(ix,2),U(ix,3),'.', ...
                    'Color',[0.87,0.43,0.18],'MarkerSize',5, ...
                    'DisplayName','Upward');
                nShown=nShown+numel(ix);
            end
            axis(ax,'equal');xlim(ax,[0,P]);ylim(ax,[0,P]);
            set(ax,'XTick',ticks,'YTick',ticks, ...
                'XTickLabel',labels,'YTickLabel',labels, ...
                'TickLabelInterpreter','tex','FontSize',12,'Layer','top');
            xlabel(ax,'\theta_2 mod 2\pi');
            ylabel(ax,'\theta_3 mod 2\pi');
            title(ax,sprintf('\\theta_1 = 0 mod 2\\pi (%s), %s', ...
                sectionDirection,tag));
            grid(ax,'on');box(ax,'on');
            if nShown==0
                text(ax,pi,pi,'No crossings in selected orientation', ...
                    'HorizontalAlignment','center');
            elseif strcmpi(sectionDirection,'both')
                legend(ax,'show','Location','northeast','AutoUpdate','off');
            end
        case 'phase'
            fig=gcf; ax=nexttile;
            nStored=size(r.cloud,2);
            if isinf(maxPhase)
                ix=1:nStored;
            else
                % IMPORTANT: a CONTIGUOUS portion, not widely spaced
                % decimation (which would make curves polygonal).
                nTake=min(nStored,max(2,floor(maxPhase)));
                ix=(nStored-nTake+1):nStored;
            end
            X=mod(r.cloud(:,ix),P);
            breaks=any(abs(diff(X,1,2))>pi,1);  % each cube face
            % Insert NaN columns AFTER each face crossing; MATLAB breaks
            % lines there instead of connecting opposite identified faces.
            n=size(X,2);
            pos=(1:n)+[0,cumsum(breaks)];
            Xline=nan(3,n+nnz(breaks));
            Xline(:,pos)=X;
            plot3(ax,Xline(1,:),Xline(2,:),Xline(3,:),'-', ...
                'LineWidth',phaseWidth,'Color',phaseColor);
            xlim(ax,[0,P]);ylim(ax,[0,P]);zlim(ax,[0,P]);
            daspect(ax,[1,1,1]);pbaspect(ax,[1,1,1]);
            set(ax,'XTick',ticks,'YTick',ticks,'ZTick',ticks, ...
                'XTickLabel',labels,'YTickLabel',labels, ...
                'ZTickLabel',labels,'TickLabelInterpreter','tex', ...
                'FontSize',11);
            xlabel(ax,'\theta_1');
            ylabel(ax,'\theta_2');
            zlabel(ax,'\theta_3');
            grid(ax,'on');box(ax,'on');view(ax,viewAzimuth,viewElevation);
            title(ax,['Orbit in the identification cube, ' tag]);
            fprintf(['  Phase plot: %d consecutive samples out of %d saved; ' ...
                '%d cube-face line breaks.\n'],n,nStored,nnz(breaks));
        otherwise
            error('Unknown plot kind %s',kind);
    end
end

function ix = fixed_plot_indices(n,maxPoints)
    if n==0, ix=[];
    elseif n<=maxPoints, ix=1:n;
    else, ix=unique(round(linspace(1,n,max(1,floor(maxPoints)))));
    end
end

function fixed_log(~,fmt,varargin)
fprintf([fmt '\n'],varargin{:});
end
