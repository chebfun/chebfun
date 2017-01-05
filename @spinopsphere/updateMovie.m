function options = updateMovie(S, dt, p, options, t, v, compGrid, plotGrid)
%UPDATEMOVIE   Update the movie when solving a PDE specified by a SPINOPSPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
nVars = S.numVars;
Clim = options{1};
dataToPlot = options{3};
ll = compGrid{1};
tt = compGrid{2};
N = size(ll, 1);
lll = plotGrid{1};
ttt = plotGrid{2};
Nplot = size(lll, 1);

for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vv = dataToPlot(v(idx:idx+N-1,:));
    
    % Change axes if necessary:
    if ( nargout == 1 )
        minvnew = min(vv(:));
        maxvnew = max(vv(:));
        if ( maxvnew > Clim(2*(k-1) + 2) )
            vscalenew = max(abs(minvnew), maxvnew);
            Clim(2*(k-1) + 2) = maxvnew + .1*vscalenew;
        end
        if ( minvnew < Clim(2*(k-1) + 1) )
            vscalenew = max(abs(minvnew), maxvnew);
            Clim(2*(k-1) + 1) = minvnew - .1*vscalenew;
        end
    end
    
    % Interpolate each variable on a finer grid:
    if ( Nplot > N )
        vvv = interp2(ll, tt, vv, lll, ttt, 'spline');
    else
        vvv = vv;
    end
    
    % Update each variable:
    vvv = vvv([floor(Nplot/2)+1:Nplot 1], :);
    set(p{k}, 'cdata', vvv)
    set(p{k}.Parent, 'clim', [Clim(2*(k-1) + 1), Clim(2*(k-1) + 2)])
    drawnow
    
end

% Update title:
titleString = sprintf('n = m = %i (DoFs = %i), dt = %1.1e, t = %.4f', N, ...
    nVars*N^2, dt, t);
set(p{nVars + 1}, 'String', titleString)

% Update outputs:
options{1} = Clim;

end