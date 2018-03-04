function [p, opts] = initializeMovie(S, dt, pref, v, compGrid, plotGrid)
%INITIALIZEMOVIE   Initialize a movie when solving a PDE specified by a SPINOP2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note: P is a NVARSx2 CELL-ARRAY that stores the NVARS plots in the first row
% and the NVARS titles in the second row. OPTS is a 3x1 CELL-ARRAY that stores
% the limits of the colorbar in OPTS{1}, the viewpoint specification in OPTS{2}
% and what kind of data to plot in OPTS{3} (real/imag/abs; when the data is
% complex-valued).

% Set-up:
dom = S.domain;                       % Spatial domain
nVars = S.numVars;                    % Number of variables (>1 for systems)
vscale = max(abs(v(:)));              % Scale of the solution
dataplot = str2func(pref.dataplot);   % Plot 'abs', 'real' or 'imag'
xx = compGrid{1};                     % Computation grid (x-direction)
yy = compGrid{2};                     % Computation grid (y-direction)
N = size(xx, 1) - 1;                  % Size of computation grid (same in x&y)
xxx = plotGrid{1};                    % Movie grid (x-direction)
yyy = plotGrid{2};                    % Movie grid (y-direction)
Nplot = size(xxx, 1) - 1;             % Size of movie grid (same in x&y)
FS = 'fontsize'; fs = 12;             % Fontsize for title

% Viewpoint specification (see SPINPREF2):
viewSpec = pref.view;
defaultPref = spinpref2();
defaultView = defaultPref.view;
while ( length(viewSpec) < 2*nVars ) 
    viewSpec = [viewSpec, defaultView];
end

% Loop over the variables:
p = cell(2, nVars); clf reset
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
    p{1,k} = surf(xxx, yyy, vvv, 'edgecolor', 'none', 'facecolor', 'interp');
    set(p{1,k}.Parent, 'clim', [Clim(2*(k-1) + 1), Clim(2*(k-1) + 2)])
    axis([dom(1) dom(2) dom(3) dom(4)])
    view(viewSpec(2*(k - 1) + 1 : 2*(k - 1) + 2))
    colorbar, colormap(pref.colormap)
    xlabel('x'), ylabel('y'), set(gca, FS, fs), box on
    
    % Plot each title:
    titleString = sprintf('Nx = Ny = %i (DoFs = %i), dt = %1.1e, t = %.4f', N, ...
        nVars*N^2, dt, 0);
    p{2,k} = title(titleString);
    drawnow
    
end

% Ask the user to press SPACE:
state = pause;
if ( strcmpi(state, 'on') == 1 )
    disp('Type <space> when ready.')
end
shg, pause

% Outputs:
opts{1} = Clim;
opts{2} = viewSpec;
opts{3} = dataplot;

end