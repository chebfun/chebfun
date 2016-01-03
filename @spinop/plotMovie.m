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

% Loop over the variables:
xxplot = [xx; 2*xx(end) - xx(end-1)];
for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vplot = dataToPlot(v(idx:idx+N-1));
    vplot = [vplot; vplot(1)]; %#ok<*AGROW>
    
    % Change axes if necessary:
    if ( nargout == 1 )
        minvnew = min(dataToPlot(vplot));
        maxvnew = max(dataToPlot(vplot));
        if ( maxvnew > Ylim(2*(k-1) + 2) )
            vscalenew = max(abs(minvnew), maxvnew);
            Ylim(2*(k-1) + 2) = maxvnew + .1*vscalenew;
        end
        if ( minvnew < Ylim(2*(k-1) + 1) )
            vscalenew = max(abs(minvnew), maxvnew);
            Ylim(2*(k-1) + 1) = minvnew - .1*vscalenew;
        end
    end
    
    % Update each variable:
    set(p{k}, 'xdata', xxplot), set(p{k}, 'ydata', vplot)
    set(p{k}.Parent, 'xlim', [dom(1), dom(2)])
    set(p{k}.Parent, 'ylim', [Ylim(2*(k-1) + 1), Ylim(2*(k-1) + 2)])
    
    % Update title:
    if ( k == 1 )
        title(p{k}.Parent, ...
            sprintf('N = %i, dt = %1.1e, t = %.4f', N, dt, t))
    end
    drawnow
    
end
plotOptions{1} = Ylim;

end