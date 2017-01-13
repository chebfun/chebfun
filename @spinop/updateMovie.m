function opts = updateMovie(S, dt, p, opts, t, v, compGrid, plotGrid)
%UPDATEMOVIE   Update the movie when solving a PDE specified by a SPINOP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
nVars = S.numVars;
Ylim = opts{1};
dataToPlot = opts{2};
xx = compGrid{1}; 
xxx = plotGrid{1};
N = size(xx, 1) - 1;
Nplot = size(xxx, 1) - 1;

for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vv = dataToPlot(v(idx:idx+N-1));
    vv = [vv; vv(1)]; %#ok<*AGROW> add repeated values (periodic endpoints)
    
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
    if ( Nplot > N )
        vvv = interp1(xx, vv, xxx, 'spline');
    else
        vvv = vv;
    end
    
    % Update each variable:
    set(p{1,k}, 'ydata', vvv)
    set(p{1,k}.Parent, 'ylim', [Ylim(2*(k-1) + 1), Ylim(2*(k-1) + 2)])
    
    % Update each title:
    titleString = sprintf('N = %i (DoFs = %i), dt = %1.1e, t = %.4f', N, ...
        nVars*N, dt, t);
    set(p{2,k}, 'String', titleString)
    drawnow
    
end

% Update outputs:
opts{1} = Ylim;

end