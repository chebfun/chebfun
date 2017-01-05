function [p, options] = initializeMovie(S, dt, pref, v, compGrid, plotGrid)
%INITIALIZEMOVIE   Initialize a movie when solving a PDE specified by a SPINOPSPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
nVars = S.numVars;
viewSpec = pref.view;
vscale = max(abs(v(:)));
dataplot = str2func(pref.dataplot);
defaultPref = spinpref2();
defaultView = defaultPref.view;
while ( length(viewSpec) < 2*nVars )
    viewSpec = [viewSpec, defaultView];
end
ll = compGrid{1};
tt = compGrid{2};
N = size(ll, 1);
lll = plotGrid{1};
ttt = plotGrid{2};
Nplot = size(lll, 1);
FS = 'fontsize';
fs = 12;

% Loop over the variables:
p = cell(nVars + 1, 1); clf reset
for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vv = dataplot(v(idx:idx+N-1,:));
    
    % Get the CLIM for the colorbar:
    if ( isempty(pref.Clim) == 1 )
        Clim(2*(k-1) + 1) = min(vv(:)) - .1*vscale; %#ok<*AGROW>
        Clim(2*(k-1) + 2) = max(vv(:)) + .1*vscale;
    else
        Clim = pref.Clim;
    end
    
    % Interpolate each variable on a finer grid:
    if ( Nplot > N )
        vvv = interp2(ll, tt, vv, lll, ttt, 'spline');
    else
        vvv = vv;
    end
    
    % Plot each variable:
    subplot(1, nVars, k)
    vvv = vvv([floor(Nplot/2)+1:Nplot 1], :);
    lam = linspace(-pi, pi, Nplot);
    th = linspace(0, pi, Nplot/2+1);
    [lam, th] = meshgrid(lam, th);
    th = pi/2 - th;
    [xxx, yyy, zzz] = sph2cart(lam, th, ones(size(lam)));
    p{k} = surf(xxx, yyy, zzz, vvv, 'edgecolor', 'none', 'facecolor', 'interp');
    set(p{k}.Parent, 'clim', [Clim(2*(k-1) + 1), Clim(2*(k-1) + 2)])
    view(viewSpec(2*(k - 1) + 1 : 2*(k - 1) + 2)), colorbar
    axis square
    set(gca, 'xtick', [-1 0 1], 'ytick', [-1 0 1], 'ztick', [-1 0 1])
    xlabel('x'), ylabel('z'), zlabel('z'), set(gca, FS, fs), box on
    drawnow
    
end

% Title:
titleString = sprintf('n = m = %i (DoFs = %i), dt = %1.1e, t = %.4f', N, ...
    nVars*N^2, dt, 0);
set(gcf, 'NextPlot', 'add');
ax = axes;
h = title(titleString);
set(ax, 'Visible', 'off', 'HandleVisibility', 'on', FS, fs);
set(h, 'Visible', 'on', 'Position', [.47 1.01 .5])

% Ask the user to press SPACE:
state = pause;
if ( strcmpi(state, 'on') == 1 )
    disp('Type <space> when ready.')
end
shg, pause

% Outputs:
p{nVars + 1} = h;
options{1} = Clim;
options{2} = viewSpec;
options{3} = dataplot;

end