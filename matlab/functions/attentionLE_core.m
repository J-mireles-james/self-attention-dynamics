% Tables 1–2, Sections 5.2–5.3. Coupled group experiments.
% Author: J.D. Mireles James
% AI use: Created and tested with AI assistance. See README for test scope.

function attentionLE_core(cfg,jobs)
% attentionLE_core  Safe long Lyapunov calculations for grouped attention.
%
% Called by the two supplied reproduction drivers. No graphics are created.
% Pure MATLAB (R2016b or later). Each job has fields:
% name, topology ('mask' or 'heads'), epsilon, A (2x2xM), B (2x2xM),
% beta (1xM), theta0 (3M x 1), and k (# leading exponents).
%
% The code is deliberately restartable and writes one result per job.
% DO NOT rename a model while recycling its checkpoint/output file.

assert(any(strcmpi(cfg.runMode,{'computeOnly','reportOnly'})), ...
       'runMode must be computeOnly or reportOnly.');
assert(exist(cfg.outputFolder,'dir')==7,'Output directory is missing.');
assert(cfg.h>0 && cfg.qrTime>0 && cfg.observationTime>0 && ...
       cfg.transientTime>=0 && cfg.checkpointTime>0);
stepsQR = exact_steps(cfg.qrTime,cfg.h,'qrTime');
nQR = exact_steps(cfg.observationTime,cfg.qrTime,'observationTime');
nBurn = exact_steps(cfg.transientTime,cfg.h,'transientTime');
chunkBurn = max(1,round(cfg.checkpointTime/cfg.h));
chunkQR = max(1,round(cfg.checkpointTime/cfg.qrTime));
historyQR = max(1,round(cfg.historyTime/cfg.qrTime));

for j=1:numel(jobs)
    job=jobs(j);
    job.theta0=mod(job.theta0(:),2*pi);
    model=build_model(job);
    assert(job.k>=1 && job.k<=model.n,'Invalid number of tangent vectors.');
    meta = struct('name',job.name,'topology',job.topology, ...
           'epsilon',job.epsilon,'A',job.A,'B',job.B, ...
           'beta',job.beta(:).','theta0',job.theta0, ...
           'h',cfg.h,'transientTime',cfg.transientTime, ...
           'observationTime',cfg.observationTime, ...
           'qrTime',cfg.qrTime,'historyTime',cfg.historyTime,'k',job.k);
    resultFile = fullfile(cfg.outputFolder,[job.name '.mat']);
    checkpointFile = fullfile(cfg.outputFolder,[job.name '_checkpoint.mat']);

    if strcmpi(cfg.runMode,'reportOnly')
        continue;
    end

    fprintf('\n=== %s: %d particles, %s, epsilon=%.6g, k=%d ===\n', ...
        job.name,model.n,job.topology,job.epsilon,job.k);
    if exist(resultFile,'file')==2 && ~cfg.forceRecompute
        tmp=load(resultFile,'meta');
        assert(isequaln(tmp.meta,meta), ...
            'Existing result %s has incompatible parameters.',resultFile);
        fprintf('Already completed. Use reportOnly for output, or forceRecompute=true.\n');
        continue;
    end
    if cfg.forceRecompute && exist(checkpointFile,'file')==2
        delete(checkpointFile);
    end
    if cfg.forceRecompute && exist(resultFile,'file')==2
        delete(resultFile);
    end
    if cfg.checkJacobian
        jacErr=jacobian_error(job.theta0,model);
        fprintf('Analytic Jacobian finite-difference relative error: %.3e\n',jacErr);
        assert(jacErr<2e-5,'Analytic Jacobian check failed.');
    end

    if cfg.resumeFromCheckpoint && exist(checkpointFile,'file')==2
        tmp=load(checkpointFile,'st','meta');
        assert(isequaln(tmp.meta,meta), ...
            'Checkpoint parameters differ. Choose a new job name or delete checkpoint.');
        st=tmp.st;
        fprintf('Resuming saved %s stage.\n',st.stage);
    else
        st=struct('stage','burn','x',job.theta0,'burnSteps',0, ...
            'block',0,'Q',[],'growth',[],'divSum',0, ...
            'times',zeros(1,0),'history',zeros(job.k,0), ...
            'spreadSum',zeros(1,model.M),'spreadCount',0, ...
            'halfGrowth',zeros(job.k,1));
    end

    % Burn-in: checkpoint even during the initial transient.
    if strcmp(st.stage,'burn')
        while st.burnSteps<nBurn
            nt=min(chunkBurn,nBurn-st.burnSteps);
            for q=1:nt
                st.x = rk4_state(st.x,cfg.h,model);
            end
            st.x=mod(st.x,2*pi);
            st.burnSteps=st.burnSteps+nt;
            save_checkpoint(checkpointFile,st,meta);
            fprintf('  transient %g / %g\n',st.burnSteps*cfg.h,cfg.transientTime);
        end
        st.stage='lyap';
        st.Q=eye(model.n,job.k);  % same thin-QR initialization as source Python
        st.growth=zeros(job.k,1);
        st.divSum=0;
        st.block=0;
        st.times=zeros(1,0);
        st.history=zeros(job.k,0);
        st.spreadSum=zeros(1,model.M);
        st.spreadCount=0;
        st.halfGrowth=zeros(job.k,1);
        save_checkpoint(checkpointFile,st,meta);
    end

    % Integrate state and variational equation in blocks, QR at unit time.
    for b=st.block+1:nQR
        for q=1:stepsQR
            [st.x,st.Q,inc] = rk4_variational(st.x,st.Q,cfg.h,model);
            st.divSum=st.divSum+inc;
        end
        [st.Q,R] = qr(st.Q,0);
        st.growth=st.growth+log(max(abs(diag(R)),realmin));
        st.x=mod(st.x,2*pi);
        st.block=b;
        if b==floor(nQR/2)
            st.halfGrowth=st.growth;
        end
        if b==1 || mod(b,historyQR)==0 || b==nQR
            t=b*cfg.qrTime;
            st.times(end+1)=t; %#ok<AGROW>
            st.history(:,end+1)=sort(st.growth/t,'descend'); %#ok<AGROW>
        end
        if b>floor(nQR/2) && (mod(b,historyQR)==0 || b==nQR)
            st.spreadSum=st.spreadSum+group_spread(st.x,model.M);
            st.spreadCount=st.spreadCount+1;
        end
        if mod(b,chunkQR)==0 || b==nQR
            save_checkpoint(checkpointFile,st,meta);
            le=sort(st.growth/(b*cfg.qrTime),'descend');
            fprintf('  Lyapunov t=%g / %g: lambda1=%+.7f, lambda2=%+.7f, count(>0.01)=%d\n', ...
                b*cfg.qrTime,cfg.observationTime,le(1),le(min(2,numel(le))),sum(le>0.01));
        end
    end

    vals=sort(st.growth/(nQR*cfg.qrTime),'descend');
    R=struct();
    R.name=job.name; R.topology=job.topology; R.epsilon=job.epsilon;
    R.particles=model.n; R.groups=model.M; R.k=job.k;
    R.lambda=vals(:).';
    R.countAbove001=sum(vals>0.01);
    R.secondHalfLambda=sort((st.growth-st.halfGrowth)/ ...
        ((nQR-floor(nQR/2))*cfg.qrTime),'descend').';
    R.secondHalfCountAbove001=sum(R.secondHalfLambda>0.01);
    R.sumPositive=sum(vals(vals>0));
    R.closestToZero=min(abs(vals));
    R.times=st.times;
    R.history=st.history;
    R.meanDivergence=st.divSum/cfg.observationTime;
    R.fullSpectrum=(job.k==model.n);
    if R.fullSpectrum
        R.divergenceDiscrepancy=sum(vals)-R.meanDivergence;
    else
        R.divergenceDiscrepancy=NaN;
    end
    R.meanGroupSpread=st.spreadSum/max(1,st.spreadCount);
    R.lastTheta=st.x;
    R.initialTheta=job.theta0;
    R.h=cfg.h; R.burn=cfg.transientTime; R.obs=cfg.observationTime;
    tempResult=[resultFile '.new.mat'];
    save(tempResult,'R','meta','-v7');
    [ok,msg]=movefile(tempResult,resultFile,'f');
    assert(ok,'Could not finalize result file: %s',msg);
    if exist(checkpointFile,'file')==2, delete(checkpointFile); end
    fprintf('Saved completed experiment: %s\n',resultFile);
end

write_reports(cfg,jobs);
end

function save_checkpoint(checkpointFile,st,meta)
    tmp=[checkpointFile '.new.mat'];
    save(tmp,'st','meta','-v7');
    [ok,msg]=movefile(tmp,checkpointFile,'f');
    assert(ok,'Could not finalize checkpoint: %s',msg);
end

function n=exact_steps(t,h,name)
    n=round(t/h);
    assert(abs(n*h-t)<1e-8*max(1,t), ...
        '%s must be an integer multiple of the time interval.',name);
end

function m=build_model(job)
    ng=numel(job.beta); n=3*ng;
    assert(size(job.A,1)==2 && size(job.A,2)==2 && size(job.A,3)==ng && ...
           size(job.B,1)==2 && size(job.B,2)==2 && size(job.B,3)==ng, ...
           'A and B must be 2x2xM.');
    assert(numel(job.theta0)==n,'Initial state must have 3M entries.');
    assert(job.epsilon>=0 && isfinite(job.epsilon), 'Invalid epsilon.');
    assert(any(strcmpi(job.topology,{'mask','heads'})), 'Unknown topology.');
    m.n=n; m.M=ng;
    m.topology=lower(job.topology);
    m.eps=job.epsilon;
    m.beta=repelem(job.beta(:),3);
    m.a11=repelem(reshape(job.A(1,1,:),[],1),3);
    m.a12=repelem(reshape(job.A(1,2,:),[],1),3);
    m.a21=repelem(reshape(job.A(2,1,:),[],1),3);
    m.a22=repelem(reshape(job.A(2,2,:),[],1),3);
    m.b11=repelem(reshape(job.B(1,1,:),[],1),3);
    m.b12=repelem(reshape(job.B(1,2,:),[],1),3);
    m.b21=repelem(reshape(job.B(2,1,:),[],1),3);
    m.b22=repelem(reshape(job.B(2,2,:),[],1),3);
    gid=floor((0:n-1)/3);
    m.same=bsxfun(@eq,gid(:),gid(:).');
    m.other=~m.same;
    m.maskWeight=double(m.same)+job.epsilon*double(m.other);
end

function [f,J]=attention_rhs(theta,m)
% Vectorized receiving-group-specific attention and exact Jacobian.
% The i-th row uses parameters A_{g(i)}, B_{g(i)}, beta_{g(i)}.
    c=cos(theta(:)).'; s=sin(theta(:)).';
    tc=-s; ts=c;
    Ac=bsxfun(@times,m.a11,c)+bsxfun(@times,m.a12,s);
    As=bsxfun(@times,m.a21,c)+bsxfun(@times,m.a22,s);
    Bc=bsxfun(@times,m.b11,c)+bsxfun(@times,m.b12,s);
    Bs=bsxfun(@times,m.b21,c)+bsxfun(@times,m.b22,s);
    value=bsxfun(@times,tc.',Ac)+bsxfun(@times,ts.',As);
    score=bsxfun(@times,m.beta, ...
       bsxfun(@times,c.',Bc)+bsxfun(@times,s.',Bs));
    weight=exp(bsxfun(@minus,score,max(score,[],2)));
    if strcmp(m.topology,'mask')
        P=weight.*m.maskWeight;
        P=bsxfun(@rdivide,P,sum(P,2));
        f=sum(P.*value,2);
    else
        Pin=weight.*m.same;
        Pin=bsxfun(@rdivide,Pin,sum(Pin,2));
        if m.M>1
            Pout=weight.*m.other;
            Pout=bsxfun(@rdivide,Pout,sum(Pout,2));
        else
            Pout=zeros(m.n);
        end
        f=sum(Pin.*value,2)+m.eps*sum(Pout.*value,2);
    end
    if nargout==1, return; end

    % Derivatives for sender j and receiver i.
    Atc=bsxfun(@times,m.a11,tc)+bsxfun(@times,m.a12,ts);
    Ats=bsxfun(@times,m.a21,tc)+bsxfun(@times,m.a22,ts);
    Btc=bsxfun(@times,m.b11,tc)+bsxfun(@times,m.b12,ts);
    Bts=bsxfun(@times,m.b21,tc)+bsxfun(@times,m.b22,ts);
    diValue=-(bsxfun(@times,c.',Ac)+bsxfun(@times,s.',As));
    djValue=bsxfun(@times,tc.',Atc)+bsxfun(@times,ts.',Ats);
    diScore=bsxfun(@times,m.beta, ...
       bsxfun(@times,tc.',Bc)+bsxfun(@times,ts.',Bs));
    djScore=bsxfun(@times,m.beta, ...
       bsxfun(@times,c.',Btc)+bsxfun(@times,s.',Bts));
    if strcmp(m.topology,'mask')
        [~,J]=head_jac(P,value,diValue,djValue,diScore,djScore);
    else
        [~,J]=head_jac(Pin,value,diValue,djValue,diScore,djScore);
        if m.M>1 && m.eps~=0
            [~,Jout]=head_jac(Pout,value,diValue,djValue,diScore,djScore);
            J=J+m.eps*Jout;
        end
    end
end

function [f,J]=head_jac(P,g,diG,djG,diS,djS)
    f=sum(P.*g,2);
    gf=bsxfun(@minus,g,f);
    J=P.*(djG+gf.*djS);
    J=J+diag(sum(P.*(diG+gf.*diS),2));
end

function x2=rk4_state(x,h,m)
    f1=attention_rhs(x,m);
    f2=attention_rhs(x+.5*h*f1,m);
    f3=attention_rhs(x+.5*h*f2,m);
    f4=attention_rhs(x+h*f3,m);
    x2=x+(h/6)*(f1+2*f2+2*f3+f4);
end

function [x2,Q2,dint]=rk4_variational(x,Q,h,m)
    [f1,J1]=attention_rhs(x,m); V1=J1*Q;
    [f2,J2]=attention_rhs(x+.5*h*f1,m); V2=J2*(Q+.5*h*V1);
    [f3,J3]=attention_rhs(x+.5*h*f2,m); V3=J3*(Q+.5*h*V2);
    [f4,J4]=attention_rhs(x+h*f3,m); V4=J4*(Q+h*V3);
    x2=x+(h/6)*(f1+2*f2+2*f3+f4);
    Q2=Q+(h/6)*(V1+2*V2+2*V3+V4);
    dint=(h/6)*(trace(J1)+2*trace(J2)+2*trace(J3)+trace(J4));
end

function err=jacobian_error(x,m)
    [~,J]=attention_rhs(x,m);
    d=1e-6; Jfd=zeros(m.n);
    for k=1:m.n
        dx=zeros(m.n,1); dx(k)=d;
        Jfd(:,k)=(attention_rhs(x+dx,m)-attention_rhs(x-dx,m))/(2*d);
    end
    err=norm(J-Jfd,'fro')/max(1,norm(J,'fro'));
end

function g=group_spread(x,M)
    U=reshape(exp(1i*x(:)),3,M);
    g=1-abs(mean(U,1));
end

function write_reports(cfg,jobs)
    prefix=fullfile(cfg.outputFolder,cfg.summaryPrefix);
    textFile=[prefix '.txt']; csvFile=[prefix '.csv'];
    ft=fopen(textFile,'w'); fc=fopen(csvFile,'w');
    assert(ft>0 && fc>0,'Could not open output summary files.');
    cleaner=onCleanup(@()close_files(ft,fc)); %#ok<NASGU>
    maxK=max([jobs.k]);
    fprintf(fc,'case,topology,epsilon,groups,particles,h,burn,observation,k,positive_gt_0p01,positive_second_half_gt_0p01,closest_to_zero,sum_positive,mean_divergence,divergence_discrepancy');
    for k=1:maxK, fprintf(fc,',lambda_%02d',k); end
    fprintf(fc,'\n');
    fprintf(ft,'COUPLED ATTENTION: TABLE-READY LYAPUNOV RESULTS\n');
    fprintf(ft,'Full spectra permit divergence cross-check; thin-QR spectra do not.\n');
    fprintf(ft,'Positive count uses cutoff lambda > 0.01 (finite-time convention).\n\n');
    for j=1:numel(jobs)
        f=fullfile(cfg.outputFolder,[jobs(j).name '.mat']);
        if exist(f,'file')~=2,continue;end
        q=load(f,'R','meta');R=q.R;
        assert(isequaln(q.meta.h,cfg.h) && ...
            isequaln(q.meta.observationTime,cfg.observationTime) && ...
            isequaln(q.meta.transientTime,cfg.transientTime), ...
            'Results for %s use different numerical settings.',jobs(j).name);
        fprintf(ft,'CASE %s\n',R.name);
        fprintf(ft,'  topology=%s eps=%.6g groups=%d N=%d k=%d\n', ...
            R.topology,R.epsilon,R.groups,R.particles,R.k);
        fprintf(ft,'  RK4_h=%.7g transient=%.1f Lyapunov_time=%.1f\n',R.h,R.burn,R.obs);
        fprintf(ft,'  n_positive_above_0.01=%d, second_half_count=%d, neutral_abs_min=%.8f, sum_positive=%.8f\n', ...
            R.countAbove001,R.secondHalfCountAbove001,R.closestToZero,R.sumPositive);
        fprintf(ft,'  mean_divergence=%+.9f, full_sum_minus_divergence=%+.3e\n', ...
            R.meanDivergence,R.divergenceDiscrepancy);
        fprintf(ft,'  exponents = [');
        fprintf(ft,' %+.8f',R.lambda);
        fprintf(ft,' ]\n  second_half_exponents = [');
        fprintf(ft,' %+.8f',R.secondHalfLambda);
        fprintf(ft,' ]\n  mean_group_spread_last_half = [');
        fprintf(ft,' %.6f',R.meanGroupSpread);
        fprintf(ft,' ]\n\n');
        fprintf(fc,'%s,%s,%.9g,%d,%d,%.9g,%.8g,%.8g,%d,%d,%d,%.10g,%.10g,%.10g,%.10g', ...
            R.name,R.topology,R.epsilon,R.groups,R.particles,R.h,R.burn,R.obs, ...
            R.k,R.countAbove001,R.secondHalfCountAbove001,R.closestToZero,R.sumPositive,R.meanDivergence,R.divergenceDiscrepancy);
        for k=1:maxK
            if k<=numel(R.lambda)
                fprintf(fc,',%.11g',R.lambda(k));
            else
                fprintf(fc,',');
            end
        end
        fprintf(fc,'\n');
    end
    fprintf('Wrote copyable results: %s and %s\n',textFile,csvFile);
end

function close_files(ft,fc)
    fclose(ft);fclose(fc);
end
