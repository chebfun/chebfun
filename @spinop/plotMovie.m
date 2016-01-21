function plotOptions = plotMovie(S, dt, p, plotOptions, t, v, gridPoints)
%PLOTMOVIE   Plot a movie when solving a PDE specified by a SPINOP.
%   PLOTMOVIE

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
dom = S.domain;
xx = gridPoints;
N = size(xx, 1);
nVars = S.numVars;
Ylim = plotOptions{1};
dataToPlot = plotOptions{2};

% Grid of the computation:
xxplot = [xx; 2*xx(end) - xx(end-1)];

% Finer grid for interploation:
Nplot = 1024;
xxxplot = trigpts(Nplot, dom);
xxxplot = [xxxplot; 2*xxxplot(end) - xxxplot(end-1)];

for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vvplot = dataToPlot(v(idx:idx+N-1));
    vvplot = [vvplot; vvplot(1)]; %#ok<*AGROW>
    
    % Change axes if necessary:
    if ( nargout == 1 )
        minvnew = min(dataToPlot(vvplot));
        maxvnew = max(dataToPlot(vvplot));
        if ( maxvnew > Ylim(2*(k-1) + 2) )
            vscalenew = max(abs(minvnew), maxvnew);
            Ylim(2*(k-1) + 2) = maxvnew + .1*vscalenew;
        end
        if ( minvnew < Ylim(2*(k-1) + 1) )
            vscalenew = max(abs(minvnew), maxvnew);
            Ylim(2*(k-1) + 1) = minvnew - .1*vscalenew;
        end
    end
    
    % Interpolate each variable on a finer grid:
    vvvplot = interp1(xxplot, vvplot, xxxplot, 'spline');
    
    % Update each variable:
    set(p{k}, 'xdata', xxxplot), set(p{k}, 'ydata', vvvplot)
    set(p{k}.Parent, 'xlim', [dom(1), dom(2)])
    set(p{k}.Parent, 'ylim', [Ylim(2*(k-1) + 1), Ylim(2*(k-1) + 2)])
    drawnow
    
end

% Update title:
titleString = sprintf('Nx = %i (DoFs = %i), dt = %1.1e, t = %.4f', N, ...
    nVars*N, dt, t);
set(p{nVars + 1}, 'String', titleString)

% Update outputs:
plotOptions{1} = Ylim;

end