function [p, opts] = initializeMovie(S, dt, pref, v, compGrid, plotGrid)
%INITIALIZEMOVIE   Initialize a movie when solving a PDE specified by a 
%SPINOPSPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note: P is a NVARSx2 CELL-ARRAY that stores the NVARS plots in the first row
% and the NVARS titles in the second row. OPTS is a 3x1 CELL-ARRAY that stores
% the limits of the colorbar in OPTS{1}, the viewpoint specification in OPTS{2}
% and what kind of data to plot in OPTS{3} (real/imag/abs; when the data is
% complex-valued).

% Set-up:
nVars = S.numVars;                    % Number of variables (>1 for systems)
vscale = max(abs(v(:)));              % Scale of the solution 
dataplot = str2func(pref.dataplot);   % Plot 'abs', 'real' or 'imag'
ll = compGrid{1};                     % Computation grid (lambda-direction)
tt = compGrid{2};                     % Computation grid (theta-direction)
N = 2*(size(ll, 1) - 1);              % Size of computation grid (same in la&th)
lll = plotGrid{1};                    % Movie grid (lambda-direction)
ttt = plotGrid{2};                    % Movie grid (theta-direction)
Nplot = 2*(size(lll, 1) - 1);         % Size of movie grid (same in la&th)
FS = 'fontsize'; fs = 12;             % Fontsize for title

% Viewpoint specification (see SPINPREFSPHERE):
viewSpec = pref.view;
defaultPref = spinprefsphere();
defaultView = defaultPref.view;
while ( length(viewSpec) < 2*nVars )
    viewSpec = [viewSpec, defaultView];
end

% Meshgrid for plotting lines of longitude and latitude and latitude:
if ( strcmpi(pref.grid, 'on') == 1 )
    gridLineType = 'k-';
    LW = 'linewidth'; lw = 1;
    numGridPtsLam = 12; % Number of longitude circles
    numGridPtsTh = 6;   % Number of latitude circles
    minPlotNum = 200;   % Number of points on each circle
    [llgl, ttgl] = meshgrid(linspace(-pi, pi, numGridPtsLam + 1), ...
        linspace(0, pi, minPlotNum));
    [llgt, ttgt] = meshgrid(linspace(-pi, pi, minPlotNum), ...
        linspace(0, pi, numGridPtsTh + 1));
end
    
% Loop over the variables:
p = cell(2, nVars); clf reset
for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vv = dataplot(v(idx:idx+N-1,:));
    vv = [vv, vv(:,1)]; %#ok<*AGROW> add repeated values (periodic endpoints)
    vv = [vv; vv(1,:)];
    vv = vv([N/2+1:N 1], :); % extract values that correspond to theta in [0,pi]
    
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
    [xxx, yyy, zzz] = sph2cart(lll, pi/2 - ttt, ones(size(lll)));
    p{1,k} = surf(xxx, yyy, zzz, vvv, 'edgecolor', 'none', 'facecolor', 'interp');
    set(p{1,k}.Parent, 'clim', [Clim(2*(k-1) + 1), Clim(2*(k-1) + 2)])
    axis equal, axis off
    view(viewSpec(2*(k - 1) + 1 : 2*(k - 1) + 2))
    colorbar, colormap(pref.colormap)
    set(gca, 'xtick', [-1 0 1], 'ytick', [-1 0 1], 'ztick', [-1 0 1])
    xlabel('x'), ylabel('y'), zlabel('z'), set(gca, FS, fs), box on
    
    % Plot lines of longitude and latitude:
    if ( strcmpi(pref.grid, 'on') == 1 )
        [xxg, yyg, zzg] = sph2cart(llgl, pi/2 - ttgl, 1 + 0*llgl);
        hold on, plot3(xxg, yyg, zzg, gridLineType, LW, lw)
        [xxg, yyg, zzg] = sph2cart(llgt, pi/2 - ttgt, 1 + 0*llgt);
        plot3(xxg', yyg', zzg', gridLineType, LW, lw)
    end
    
    % Plot each title:
    titleString = sprintf('n = m = %i (DoFs = %i), dt = %1.1e, t = %.4f', N, ...
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