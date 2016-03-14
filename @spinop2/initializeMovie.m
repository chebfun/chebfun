function [p, options] = initializeMovie(S, dt, pref, v, dataGrid, plotGrid)
%INITIALIZEMOVIE   Initialize a movie when solving a PDE specified by a SPINOP2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
dom = S.domain;
nVars = S.numVars;
viewSpec = pref.view;
vscale = max(abs(v(:)));
dataToPlot = str2func(pref.dataToPlot);
defaultPref = spinpref2();
defaultView = defaultPref.view;
while ( length(viewSpec) < 2*nVars )
    viewSpec = [viewSpec, defaultView];
end
xx = dataGrid{1};
yy = dataGrid{2};
N = size(xx, 1) - 1;
xxx = plotGrid{1};
yyy = plotGrid{2};

% Loop over the variables:
p = cell(nVars + 1, 1); clf reset
for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vv = dataToPlot(v(idx:idx+N-1,:));
    vv = [vv, vv(:,1)]; %#ok<*AGROW>
    vv = [vv; vv(1,:)];
    
    % Get the CLIM for the colorbar:
    if ( isempty(pref.Clim) == 1 )
        Clim(2*(k-1) + 1) = min(vv(:)) - .1*vscale;
        Clim(2*(k-1) + 2) = max(vv(:)) + .1*vscale;
    else
        Clim = pref.Clim;
    end
    
    % Interpolate each variable on a finer grid:
    vvv = interp2(xx, yy, vv, xxx, yyy, 'spline');
    
    % Plot each variable:
    subplot(1, nVars, k)
    p{k} = surf(xxx, yyy, vvv, 'edgecolor', 'none');
    set(p{k}.Parent, 'clim', [Clim(2*(k-1) + 1), Clim(2*(k-1) + 2)])
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
options{1} = Clim;
options{2} = viewSpec;
options{3} = dataToPlot;

end