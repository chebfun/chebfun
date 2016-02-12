function [p, plotOptions] = initializeMovie(S, dt, pref, v, gridPoints)
%INITIALIZEMOVIE   Initialize a movie when solving a PDE specified by a SPINOP2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
dom = S.domain;
nVars = S.numVars;
viewSpec = pref.view;
dataToPlot = str2func(pref.dataToPlot);
defaultPref = spinpref2();
defaultView = defaultPref.view;
while ( length(viewSpec) < 2*nVars )
    viewSpec = [viewSpec, defaultView];
end
xx = gridPoints{1};
yy = gridPoints{2};
N = size(xx, 1);

% Grid of the computation:
xxplot = [xx, 2*xx(:,end) - xx(:,end-1)];
xxplot = [xxplot; xxplot(1,:)];
yyplot = [yy; 2*yy(end,:) - yy(end-1,:)];
yyplot = [yyplot, yyplot(:,1)];

% Finer grid for interploation:
Nplot = max(N, 256);
[xxxplot, yyyplot] = meshgrid(trigpts(Nplot, dom));
xxxplot = [xxxplot, 2*xxxplot(:,end) - xxxplot(:,end-1)];
xxxplot = [xxxplot; xxxplot(1,:)];
yyyplot = [yyyplot; 2*yyyplot(end,:) - yyyplot(end-1,:)];
yyyplot = [yyyplot, yyyplot(:,1)];

% Loop over the variables:
p = cell(nVars + 1, 1); clf reset
for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vvplot = dataToPlot(v(idx:idx+N-1,:));
    vvplot = [vvplot, vvplot(:,1)]; %#ok<*AGROW>
    vvplot = [vvplot; vvplot(1,:)];
    
    % Interpolate each variable on a finer grid:
    vvvplot = interp2(xxplot, yyplot, vvplot, xxxplot, yyyplot, 'spline');
    
    % Plot each variable:
    subplot(1, nVars, k)
    p{k} = surf(xxxplot, yyyplot, vvvplot, 'edgecolor', 'none');
    axis([dom(1) dom(2) dom(3) dom(4)])
    view(viewSpec(2*(k - 1) + 1 : 2*(k - 1) + 2)), colorbar
    xlabel('x'), ylabel('y'), set(gca, 'FontSize', 16), box on
    drawnow
    
end

% Title:
titleString = sprintf('Nx = Ny = %i (DoFs = %i), dt = %1.1e, t = %.4f', N, ...
    nVars*N^2, dt, 0);
set(gcf, 'NextPlot', 'add');
ax = axes;
h = title(titleString);
set(ax, 'Visible', 'off', 'HandleVisibility', 'off', 'Fontsize', 16);
set(h, 'Visible', 'on', 'Position', [.47 1.01 .5])

% Ask the user to press SPACE:
state = pause;
if ( strcmpi(state, 'on') == 1 )
    disp('Type <space> when ready.')
end
shg, pause

% Outputs:
p{nVars + 1} = h;
plotOptions{1} = viewSpec;
plotOptions{2} = dataToPlot;

end