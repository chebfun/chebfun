function [p, options] = initializeMovie(S, dt, pref, v, dataGrid, plotGrid)
%INITIALIZEMOVIE   Initialize a movie when solving a PDE specified by a SPINOP.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
dom = S.domain;
nVars = S.numVars;
vscale = max(abs(v));
dataToPlot = str2func(pref.dataToPlot);
xx = dataGrid{1};
xxx = plotGrid{1};
N = size(xx, 1) - 1;

% Loop over the variables:
p = cell(nVars, 1); clf reset
for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vvplot = dataToPlot(v(idx:idx+N-1));
    vvplot = [vvplot; vvplot(1)]; %#ok<*AGROW>
    
    % Get the YLIM for the y-axis:
    if ( isempty(pref.Ylim) == 1 )
        Ylim(2*(k-1) + 1) = min(vvplot) - .1*vscale;
        Ylim(2*(k-1) + 2) = max(vvplot) + .1*vscale;
    else
        Ylim = pref.Ylim;
    end
    
    % Interpolate each variable on a finer grid:
    vvvplot = interp1(xx, vvplot, xxx, 'spline');
    
    % Plot each variable:
    subplot(1, nVars, k)
    p{k} = plot(xxx, vvvplot, 'linewidth', 3);
    axis([dom(1), dom(2), Ylim(2*(k-1) + 1), Ylim(2*(k-1) + 2)])
    if ( nVars == 1 )
        xlabel('x'), ylabel('u(t,x)'), grid on
    else
        xlabel('x'), ylabel(['u_',num2str(k),'(t,x)']), grid on
    end
    set(gca, 'FontSize', 16), box on
    drawnow
    
end

% Title:
titleString = sprintf('N = %i (DoFs = %i), dt = %1.1e, t = %.4f', N, ...
    nVars*N, dt, 0);
set(gcf, 'NextPlot', 'add');
ax = axes;
h = title(titleString);
set(ax, 'Visible', 'off', 'HandleVisibility', 'off', 'Fontsize', 16);
set(h, 'Visible', 'on', 'Position', [.5 1.00 .5])

% Ask the user to press SPACE:
state = pause;
if ( strcmpi(state, 'on') == 1 )
    disp('Type <space> when ready.')
end
shg, pause

% Outputs:
p{nVars + 1} = h;
options{1} = Ylim;
options{2} = dataToPlot;

end