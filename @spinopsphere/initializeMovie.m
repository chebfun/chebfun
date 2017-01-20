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
nVars = S.numVars;
viewSpec = pref.view;
vscale = max(abs(v(:)));
dataplot = str2func(pref.dataplot);
defaultPref = spinprefsphere();
defaultView = defaultPref.view;
while ( length(viewSpec) < 2*nVars )
    viewSpec = [viewSpec, defaultView];
end
ll = compGrid{1};
tt = compGrid{2};
N = 2*(size(ll, 1) - 1);
lll = plotGrid{1};
ttt = plotGrid{2};
Nplot = 2*(size(lll, 1) - 1);
FS = 'fontsize';
fs = 12;

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