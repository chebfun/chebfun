function [p, opts] = initializeMovie(S, dt, pref, v, compGrid, plotGrid)
%INITIALIZEMOVIE   Initialize a movie when solving a PDE specified by a SPINOP3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note: P is a (NVARS+1)x1 CELL-ARRAY that stores the NVARS plots in its first 
% NVARS entries, and the title in its (NVARS+1)st entry. OPTS is a 3x1
% CELL-ARRAY that stores the limits of the colorbar in OPTS{1}, the slices of
% the volumetric slice plot in OPTS{2} and what kind of data to plot in OPTS{3} 
% (real/imag/abs; when the data is complex-valued).

% Set-up:
nVars = S.numVars;
vscale = max(abs(v(:)));
dataplot = str2func(pref.dataplot);
dom = S.domain;
xx = compGrid{1};
yy = compGrid{2};
zz = compGrid{3};
N = size(xx, 1) - 1;
xxx = plotGrid{1};
yyy = plotGrid{2};
zzz = plotGrid{3};
Nplot = size(xxx, 1) - 1;
FS = 'fontsize';
fs = 12;

% Slices:
ttx = trigpts(N, dom(1:2));
tty = trigpts(N, dom(3:4));
ttz = trigpts(N, dom(5:6));
ttx = [ttx; 2*ttx(end) - ttx(end-1)];
tty = [tty; 2*tty(end) - tty(end-1)];
ttz = [ttz; 2*ttz(end) - ttz(end-1)];
if ( isempty(pref.slices) == 1 )
    Sx = ttx(floor(N/2) + 1);
    Sy = tty(floor(N/2) + 1);
    Sz = ttz(floor(N/2) + 1);
else
    slices = pref.slices;
    Sx = slices{1};
    for k = 1:length(Sx)
        pos = Sx(k);
        [~, id] = min(abs(ttx - pos));
        Sx(k) = ttx(id);
    end
    Sy = slices{2};
    for k = 1:length(Sy)
        pos = Sy(k);
        [~, id] = min(abs(tty - pos));
        Sy(k) = tty(id);
    end
    Sz = slices{3};
    for k = 1:length(Sz)
        pos = Sz(k);
        [~, id] = min(abs(ttz - pos));
        Sz(k) = ttz(id);
    end
end

% Loop over the variables:
p = cell(nVars + 1, 1); clf reset
for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vv = dataplot(v(idx:idx+N-1,:,:));
    vv = [vv, vv(:,1,:)]; %#ok<*AGROW> add repeated values (periodic endpoints)
    vv = [vv; vv(1,:,:)];
    vv = cat(3, vv, vv(:,:,1));
    
    % Get the CLIM for the colorbar:
    if ( isempty(pref.Clim) == 1 )
        Clim(2*(k-1) + 1) = min(vv(:)) - .1*vscale;
        Clim(2*(k-1) + 2) = max(vv(:)) + .1*vscale;
    else
        Clim = pref.Clim;
    end
    
    % Interpolate each variable on a finer grid:
    if ( Nplot > N )
        vvv = interp3(xx, yy, zz, vv, xxx, yyy, zzz, 'spline');
    else
        vvv = vv;
    end
    
    % Plot each variable:
    subplot(1, nVars, k) 
    p{k} = slice(xxx, yyy, zzz, vvv, Sx, Sy, Sz);
    set(p{k}, 'edgecolor', 'none', 'facecolor', 'interp')
    ax = p{k}.Parent; set(ax, 'clim', [Clim(2*(k-1) + 1), Clim(2*(k-1) + 2)])
    axis equal, axis([dom(1) dom(2) dom(3) dom(4) dom(5) dom(6)])
    colorbar, colormap(pref.colormap)
    xlabel('x'), ylabel('y'), zlabel('z'), set(gca, FS, fs), box on
    drawnow
    
end

% Title:
titleString = sprintf('Nx = Ny = Nz = %i (DoFs = %i), dt = %1.1e, t = %.4f', ...
    N, nVars*N^3, dt, 0);
set(gcf, 'NextPlot', 'add');
ax = axes;
h = title(titleString);
set(ax, 'Visible', 'off', 'HandleVisibility', 'on', FS, fs)
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
opts{2} = {Sx, Sy, Sz};
opts{3} = dataplot;

end