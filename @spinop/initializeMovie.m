function [p, opts] = initializeMovie(S, dt, pref, v, compGrid, plotGrid)
%INITIALIZEMOVIE   Initialize a movie when solving a PDE specified by a SPINOP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note: P is a (NVARS+1)x1 CELL-ARRAY that stores the NVARS plots in its first 
% NVARS entries, and the title in its (NVARS+1)st entry. OPTS is a 2x1
% CELL-ARRAY that stores the limits of the y-axis in OPTS{1} and what kind of
% data to plot in OPTS{2} (real/imag/abs; when the data is complex-valued).

% Set-up:
dom = S.domain;
nVars = S.numVars;
vscale = max(abs(v));
dataplot = str2func(pref.dataplot);
xx = compGrid{1};
xxx = plotGrid{1};
N = size(xx, 1) - 1;
Nplot = size(xxx, 1) - 1;
FS = 'fontsize';
fs = 12;

% Loop over the variables:
p = cell(nVars, 1); clf reset
for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vv = dataplot(v(idx:idx+N-1));
    vv = [vv; vv(1)]; %#ok<*AGROW>
    
    % Get the YLIM for the y-axis:
    if ( isempty(pref.Ylim) == 1 )
        Ylim(2*(k-1) + 1) = min(vv) - .1*vscale;
        Ylim(2*(k-1) + 2) = max(vv) + .1*vscale;
    else
        Ylim = pref.Ylim;
    end
    
    % Interpolate each variable on a finer grid:
    if ( Nplot > N )
        vvv = interp1(xx, vv, xxx, 'spline');
    else
        vvv = vv;
    end
    
    % Plot each variable:
    subplot(1, nVars, k)
    p{k} = plot(xxx, vvv, 'linewidth', 3);
    axis([dom(1), dom(2), Ylim(2*(k-1) + 1), Ylim(2*(k-1) + 2)])
    if ( nVars == 1 )
        xlabel('x'), ylabel('u(t,x)'), grid on
    else
        xlabel('x'), ylabel(['u_',num2str(k),'(t,x)']), grid on
    end
    set(gca, FS, fs), box on
    drawnow
    
end

% Title:
titleString = sprintf('N = %i (DoFs = %i), dt = %1.1e, t = %.4f', N, ...
    nVars*N, dt, 0);
set(gcf, 'NextPlot', 'add');
ax = axes;
h = title(titleString);
set(ax, 'Visible', 'off', 'HandleVisibility', 'on', FS, fs);
set(h, 'Visible', 'on', 'Position', [.5 1.00 .5])

% Ask the user to press SPACE:
state = pause;
if ( strcmpi(state, 'on') == 1 )
    disp('Type <space> when ready.')
end
shg, pause

% Outputs:
p{nVars + 1} = h;
opts{1} = Ylim;
opts{2} = dataplot;

end