% Sections 3.1 and 3.3. Periodic orbits in the identification cube.
% Author: J.D. Mireles James
% AI use: Created and tested with AI assistance. See README for test scope.

function portrait(data,recompute)
c = data.caseSpec;
if recompute
    rhs = @(~,x) portrait_field(x,c.A,c.B,c.beta);
    opts = odeset('RelTol',data.relativeTolerance,'AbsTol',data.absoluteTolerance,'MaxStep',data.maxStep);
    [~,y] = ode113(rhs,data.tOrbit,double(c.thetaSeed(:)),opts);
    data.periodicLift = y.';
    for j = 1:size(data.transientStarts,2)
        [~,y] = ode113(rhs,data.tTransient,data.transientStarts(:,j),opts);
        data.transientLift{j} = y.';
    end
end
figure; tiledlayout(1,2);
for show = [false true]
    portrait_plot(data,2,0.6,Inf,true,35,25,show,true);
end
end
function f = portrait_field(theta,A,B,beta)
% Same autonomous 3-particle softmax vector field as the bifurcation code.
    theta = theta(:);
    u = [cos(theta).';sin(theta).'];
    v = [-u(2,:);u(1,:)];
    scores = beta*(u.'*B*u);
    scores = scores-max(scores,[],2); % stable row-wise softmax
    weights = exp(scores);
    weights = weights./sum(weights,2);
    avg = u*weights.';
    f = sum(v.*(A*avg),1).';
end

function portrait_plot(data,orbitWidth,transWidth,maxTransient, ...
        showCube,az,el,showTransients,showStartMarkers)
% Graphics intentionally simple: ~10 line objects, no surface meshes.
    c = data.caseSpec;
    ax = nexttile;
    hold(ax,'on');

    if showCube
        portrait_cube_edges(ax);
    end

    if showTransients
        n = numel(data.transientLift);
        transColor = [0.34 0.37 0.41];
        for j = 1:n
            X = data.transientLift{j};
            % Keep a contiguous block. Never connect distant samples.
            if isfinite(maxTransient) && size(X,2)>maxTransient
                X = X(:,1:floor(maxTransient));
            end
            Z = portrait_wrap_with_nan(X);
            plot3(ax,Z(1,:),Z(2,:),Z(3,:),'-', ...
                'Color',transColor,'LineWidth',transWidth, ...
                'HandleVisibility','off');
        end
        if showStartMarkers && ~isempty(data.transientStarts)
            S = mod(data.transientStarts,2*pi);
            plot3(ax,S(1,:),S(2,:),S(3,:),'o', ...
                'LineStyle','none','Color',transColor, ...
                'MarkerFaceColor',transColor,'MarkerSize',4, ...
                'HandleVisibility','off');
        end
    end

    % Plot the periodic orbit LAST so it remains visually prominent.
    Z = portrait_wrap_with_nan(data.periodicLift);
    plot3(ax,Z(1,:),Z(2,:),Z(3,:),'-', ...
        'Color',c.orbitColor,'LineWidth',orbitWidth);

    L=2*pi;
    xlim(ax,[0,L]); ylim(ax,[0,L]); zlim(ax,[0,L]);
    daspect(ax,[1,1,1]);
    view(ax,az,el);
    ticks=[0,pi,2*pi]; labels={'0','\pi','2\pi'};
    set(ax,'XTick',ticks,'YTick',ticks,'ZTick',ticks, ...
        'XTickLabel',labels,'YTickLabel',labels,'ZTickLabel',labels, ...
        'FontSize',11);
    xlabel(ax,'\theta_1'); ylabel(ax,'\theta_2'); zlabel(ax,'\theta_3');
    title(ax,c.titleText,'Interpreter','tex');
    grid(ax,'on'); box(ax,'on');
end

function Z = portrait_wrap_with_nan(X)
% Insert NaNs between points straddling identified faces of [0,2pi]^3.
    if isempty(X), Z = X; return; end
    Y = mod(X,2*pi);
    if size(Y,2)<2, Z=Y; return; end
    jumps = any(abs(diff(Y,1,2))>pi,1);
    n = size(Y,2);
    position = (1:n) + [0,cumsum(jumps)];
    Z = nan(3,n+sum(jumps));
    Z(:,position) = Y;
end

function portrait_cube_edges(ax)
% One simple wireframe cube, no transparent surfaces or lighting.
    L = 2*pi;
    V = L*[0 0 0; 1 0 0; 1 1 0; 0 1 0; ...
           0 0 1; 1 0 1; 1 1 1; 0 1 1].';
    E = [1 2;2 3;3 4;4 1;5 6;6 7;7 8;8 5; ...
         1 5;2 6;3 7;4 8];
    coords = nan(3,3*size(E,1));
    for j=1:size(E,1)
        coords(:,3*j-2)=V(:,E(j,1));
        coords(:,3*j-1)=V(:,E(j,2));
    end
    plot3(ax,coords(1,:),coords(2,:),coords(3,:),'-', ...
        'Color',[0.6 0.6 0.6],'LineWidth',0.8, ...
        'HandleVisibility','off');
end
