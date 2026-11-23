% Figure 9, Section 5.1. Six-dimensional RK4 and QR accumulation.
% Author: J.D. Mireles James
% AI use: Created and tested with AI assistance. See README for test scope.

function R = coupled_compute_case(cfg,ep,kind,nBurn,nBlocks,stepsQr, ...
                      historyBlocks,checkpointBlocks,ckptFile,doResume)
    if doResume && exist(ckptFile,'file')
        S = load(ckptFile,'st','configSaved','savedEps','savedKind');
        coupled_check_config(S.configSaved,cfg);
        assert(S.savedEps==ep && strcmp(S.savedKind,kind), ...
            'Checkpoint model differs from selected experiment.');
        st=S.st;
        fprintf('  Resuming QR block %d / %d.\n',st.b,nBlocks);
    else
        tic;
        fprintf('  Removing initial transient: %g time units...\n', ...
                cfg.transientTime);
        theta = cfg.initialTheta;
        for k=1:nBurn
            theta = coupled_rk4_state(theta,cfg.h,cfg.A,cfg.B, ...
                                      cfg.beta,ep,kind);
            if mod(k,round(100/cfg.h))==0
                theta=mod(theta,2*pi); % same point on T^6; prevents angle growth
            end
            if mod(k,round(250/cfg.h))==0
                fprintf('  Transient t=%g/%g, elapsed %.1fs\n', ...
                    k*cfg.h,cfg.transientTime,toc);
            end
        end
        st=struct('theta',mod(theta,2*pi),'Q',eye(6), ...
            'growth',zeros(6,1),'divIntegral',0,'b',0, ...
            'times',zeros(1,0),'history',zeros(6,0));
        fprintf('  Transient removed in %.1f seconds.\n',toc);
    end
    tic;
    for b=st.b+1:nBlocks
        for k=1:stepsQr
            [st.theta,st.Q,d] = coupled_rk4_variational( ...
                st.theta,st.Q,cfg.h,cfg.A,cfg.B,cfg.beta,ep,kind);
            st.divIntegral=st.divIntegral+d;
        end
        [st.Q,Rmat]=qr(st.Q,0);
        st.growth=st.growth+log(max(abs(diag(Rmat)),realmin));
        st.theta=mod(st.theta,2*pi);
        st.b=b;
        if b==1 || mod(b,historyBlocks)==0 || b==nBlocks
            tt=b*cfg.qrTime;
            le=sort(st.growth/tt,'descend');
            st.times(end+1)=tt; %#ok<AGROW>
            st.history(:,end+1)=le; %#ok<AGROW>
        end
        if mod(b,checkpointBlocks)==0 || b==nBlocks
            savedEps=ep; savedKind=kind; configSaved=cfg; %#ok<NASGU>
            % Progress is reported in memory.
            currentLE=sort(st.growth/(b*cfg.qrTime),'descend');
            fprintf('  t_LE=%6g / %g, lambda1=%+.5f, lambda2=%+.5f, elapsed %.1fs\n', ...
                b*cfg.qrTime,cfg.observationTime, ...
                currentLE(1),currentLE(2),toc);
        end
    end
    le=sort(st.growth/cfg.observationTime,'descend');
    R=struct('kind',kind,'epsilon',ep,'exponents',le(:).', ...
        'times',st.times,'history',st.history, ...
        'meanDivergence',st.divIntegral/cfg.observationTime, ...
        'divergenceDiscrepancy',sum(le)-st.divIntegral/cfg.observationTime, ...
        'terminalTheta',st.theta,'initialTheta',cfg.initialTheta, ...
        'transientTime',cfg.transientTime, ...
        'observationTime',cfg.observationTime);
end

function theta2 = coupled_rk4_state(theta,h,A,B,beta,ep,kind)
    k1=coupled_field_jac(theta,A,B,beta,ep,kind);
    k2=coupled_field_jac(theta+.5*h*k1,A,B,beta,ep,kind);
    k3=coupled_field_jac(theta+.5*h*k2,A,B,beta,ep,kind);
    k4=coupled_field_jac(theta+h*k3,A,B,beta,ep,kind);
    theta2=theta+(h/6)*(k1+2*k2+2*k3+k4);
end

function [theta2,Q2,divStep] = coupled_rk4_variational(theta,Q,h,A,B,beta,ep,kind)
    [f1,J1]=coupled_field_jac(theta,A,B,beta,ep,kind);
    k1=J1*Q;
    [f2,J2]=coupled_field_jac(theta+.5*h*f1,A,B,beta,ep,kind);
    k2=J2*(Q+.5*h*k1);
    [f3,J3]=coupled_field_jac(theta+.5*h*f2,A,B,beta,ep,kind);
    k3=J3*(Q+.5*h*k2);
    [f4,J4]=coupled_field_jac(theta+h*f3,A,B,beta,ep,kind);
    k4=J4*(Q+h*k3);
    theta2=theta+(h/6)*(f1+2*f2+2*f3+f4);
    Q2=Q+(h/6)*(k1+2*k2+2*k3+k4);
    divStep=(h/6)*(trace(J1)+2*trace(J2)+2*trace(J3)+trace(J4));
end

function [f,J] = coupled_field_jac(theta,A,B,beta,ep,kind)
% Analytic vector field and Jacobian on T^6. Particle index i belongs to
% group 1 for i=1:3 and group 2 for i=4:6.
% s_ij=beta u_i^T B u_j, g_ij=t_i^T A u_j.
% For a normalized head, f_i=sum_j p_ij g_ij and
% J_ij=p_ij*(d_k g_ij+(g_ij-f_i)*d_k s_ij) +
%      delta_ij*sum_k p_ik*(d_i g_ik+(g_ik-f_i)*d_i s_ik).
    c=cos(theta(:).'); s=sin(theta(:).');
    U=[c;s]; V=[-s;c];
    scores=beta*(U.'*B*U);
    g=V.'*A*U;
    same=logical([ones(3),zeros(3);zeros(3),ones(3)]);
    if strcmp(kind,'mask')
        logits=scores;
        logits(~same)=logits(~same)+log(ep);
        P=coupled_softmax(logits);
        if nargout==1
            f=sum(P.*g,2); return;
        end
    elseif strcmp(kind,'heads')
        logitsIn=scores;
        logitsIn(~same)=-Inf;
        P=coupled_softmax(logitsIn);
        logitsOut=scores;
        logitsOut(same)=-Inf;
        Pc=coupled_softmax(logitsOut);
        if nargout==1
            f=sum(P.*g,2)+ep*sum(Pc.*g,2); return;
        end
    else
        error('Unsupported coupling type: %s',kind);
    end
    ds_i=beta*(V.'*B*U);  % i angular derivative of s_ij
    dg_i=-(U.'*A*U);    % i angular derivative of g_ij
    ds_j=beta*(U.'*B*V);  % j angular derivative of s_ij
    dg_j=V.'*A*V;       % j angular derivative of g_ij
    [f,J]=coupled_head_jac(P,g,ds_i,dg_i,ds_j,dg_j);
    if strcmp(kind,'heads')
        [fc,Jc]=coupled_head_jac(Pc,g,ds_i,dg_i,ds_j,dg_j);
        f=f+ep*fc;
        J=J+ep*Jc;
    end
end

function P=coupled_softmax(logits)
    z=exp(bsxfun(@minus,logits,max(logits,[],2)));
    P=bsxfun(@rdivide,z,sum(z,2));
end

function [f,J]=coupled_head_jac(P,g,ds_i,dg_i,ds_j,dg_j)
    f=sum(P.*g,2);
    gf=bsxfun(@minus,g,f);
    self=sum(P.*(dg_i+gf.*ds_i),2);
    J=P.*(dg_j+gf.*ds_j)+diag(self);
end

function err=coupled_jacobian_test(x,A,B,beta,eps,kind)
    [~,J]=coupled_field_jac(x,A,B,beta,eps,kind);
    Jfd=zeros(6);
    delta=1e-6;
    for k=1:6
        e=zeros(6,1);e(k)=delta;
        fp=coupled_field_jac(x+e,A,B,beta,eps,kind);
        fm=coupled_field_jac(x-e,A,B,beta,eps,kind);
        Jfd(:,k)=(fp-fm)/(2*delta);
    end
    err=norm(J-Jfd,'fro')/max(1,norm(J,'fro'));
end

function coupled_check_config(c1,c2)
% Prevent loading cached data from an unrelated numerical experiment.
    keys={'A','B','beta','h','transientTime','observationTime', ...
          'qrTime','historyTime','initialTheta','epsilonList'};
    for k=1:numel(keys)
        key=keys{k};
        assert(isequaln(c1.(key),c2.(key)), ...
            'Saved data and current %s differ. Change outputPrefix or settings.',key);
    end
end
