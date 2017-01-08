function [p, opts] = initializeMovie(S, dt, pref, v, compGrid, plotGrid)
%INITIALIZEMOVIE   Initialize a movie when solving a PDE specified by a SPINOP2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note: P is a (NVARS+1)x1 CELL-ARRAY that stores the NVARS plots in its first 
% NVARS entries, and the title in its (NVARS+1)st entry. OPTS is a 3x1
% CELL-ARRAY that stores the limits of the colorbar in OPTS{1}, the viewpoint 
% specification in OPTS{2} and what kind of data to plot in OPTS{3} 
% (real/imag/abs; when the data is complex-valued).

% Set-up:
dom = S.domain;
nVars = S.numVars;
viewSpec = pref.view;
vscale = max(abs(v(:)));
dataplot = str2func(pref.dataplot);
defaultPref = spinpref2();
defaultView = defaultPref.view;
while ( length(viewSpec) < 2*nVars )
    viewSpec = [viewSpec, defaultView];
end
xx = compGrid{1};
yy = compGrid{2};
N = size(xx, 1) - 1;
xxx = plotGrid{1};
yyy = plotGrid{2};
Nplot = size(xxx, 1) - 1;
FS = 'fontsize';
fs = 12;

% Loop over the variables:
p = cell(nVars + 1, 1); clf reset
for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vv = dataplot(v(idx:idx+N-1,:));
    vv = [vv, vv(:,1)]; %#ok<*AGROW> add repeated values (periodic endpoints)
    vv = [vv; vv(1,:)]; 
    
    % Get the CLIM for the colorbar:
    if ( isempty(pref.Clim) == 1 )
        Clim(2*(k-1) + 1) = min(vv(:)) - .1*vscale;
        Clim(2*(k-1) + 2) = max(vv(:)) + .1*vscale;
    else
        Clim = pref.Clim;
    end
    
    % Interpolate each variable on a finer grid:
    if ( Nplot > N ) 
        vvv = interp2(xx, yy, vv, xxx, yyy, 'spline');
    else
        vvv = vv;
    end
    
    % Plot each variable:
    subplot(1, nVars, k)
    p{k} = surf(xxx, yyy, vvv, 'edgecolor', 'none', 'facecolor', 'interp');
    set(p{k}.Parent, 'clim', [Clim(2*(k-1) + 1), Clim(2*(k-1) + 2)])
    axis([dom(1) dom(2) dom(3) dom(4)])
    view(viewSpec(2*(k - 1) + 1 : 2*(k - 1) + 2)), colorbar
    xlabel('x'), ylabel('y'), set(gca, FS, fs), box on
    drawnow
    
end

% Title:
titleString = sprintf('Nx = Ny = %i (DoFs = %i), dt = %1.1e, t = %.4f', N, ...
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
opts{1} = Clim;
opts{2} = viewSpec;
opts{3} = dataplot;

end