function options = plotMovie(S, dt, p, options, t, v, dataGrid, plotGrid)
%PLOTMOVIE   Plot a movie when solving a PDE specified by a SPINOP.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
dom = S.domain;
nVars = S.numVars;
Ylim = options{1};
dataToPlot = options{2};
xx = dataGrid{1};
xxx = plotGrid{1};
N = size(xx, 1) - 1;

for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vv = dataToPlot(v(idx:idx+N-1));
    vv = [vv; vv(1)]; %#ok<*AGROW>
    
    % Change axes if necessary:
    if ( nargout == 1 )
        minvnew = min(dataToPlot(vv));
        maxvnew = max(dataToPlot(vv));
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
    vvv = interp1(xx, vv, xxx, 'spline');
    
    % Update each variable:
    set(p{k}, 'xdata', xxx), set(p{k}, 'ydata', vvv)
    set(p{k}.Parent, 'xlim', [dom(1), dom(2)])
    set(p{k}.Parent, 'ylim', [Ylim(2*(k-1) + 1), Ylim(2*(k-1) + 2)])
    drawnow
    
end

% Update title:
titleString = sprintf('N = %i (DoFs = %i), dt = %1.1e, t = %.4f', N, ...
    nVars*N, dt, t);
set(p{nVars + 1}, 'String', titleString)

% Update outputs:
options{1} = Ylim;

end