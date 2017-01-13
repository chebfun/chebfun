function opts = updateMovie(S, dt, p, opts, t, v, compGrid, plotGrid)
%UPDATEMOVIE   Update the movie when solving a PDE specified by a SPINOPSPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
nVars = S.numVars;
Clim = opts{1};
dataToPlot = opts{3};
ll = compGrid{1};
tt = compGrid{2};
N = 2*(size(ll, 1) - 1);
lll = plotGrid{1};
ttt = plotGrid{2};
Nplot = 2*(size(lll, 1) - 1);

for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vv = dataToPlot(v(idx:idx+N-1,:));
    vv = [vv, vv(:,1)]; %#ok<*AGROW> add repeated values (periodic endpoints)
    vv = [vv; vv(1,:)];
    vv = vv([N/2+1:N 1], :); % extract values that correspond to theta in [0,pi]
    
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
    set(p{1,k}, 'cdata', vvv)
    set(p{1,k}.Parent, 'clim', [Clim(2*(k-1) + 1), Clim(2*(k-1) + 2)])
    
    % Update title:
    titleString = sprintf('n = m = %i (DoFs = %i), dt = %1.1e, t = %.4f', N, ...
        nVars*N^2, dt, t);
    set(p{2,k}, 'String', titleString)
    drawnow

end

% Update outputs:
opts{1} = Clim;

end