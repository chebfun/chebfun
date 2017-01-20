function opts = updateMovie(S, dt, p, opts, t, v, compGrid, plotGrid)
%UPDATEMOVIE   Update the movie when solving a PDE specified by a SPINOP2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
nVars = S.numVars;
Clim = opts{1};
dataToPlot = opts{3};
xx = compGrid{1};
yy = compGrid{2};
N = size(xx, 1) - 1;
xxx = plotGrid{1};
yyy = plotGrid{2};
Nplot = size(xxx, 1) - 1;

for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vv = dataToPlot(v(idx:idx+N-1,:));
    vv = [vv, vv(:,1)]; %#ok<*AGROW> add repeated values (periodic endpoints)
    vv = [vv; vv(1,:)];
    
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
        vvv = interp2(xx, yy, vv, xxx, yyy, 'spline');
    else
        vvv = vv;
    end
    
    % Update each variable:
    set(p{1,k}, 'zdata', vvv)
    set(p{1,k}.Parent, 'clim', [Clim(2*(k-1) + 1), Clim(2*(k-1) + 2)])
    
    % Update each title:
    titleString = sprintf('Nx = Ny = %i (DoFs = %i), dt = %1.1e, t = %.4f', N, ...
        nVars*N^2, dt, t);
    set(p{2,k}, 'String', titleString)
    drawnow

end

% Update outputs:
opts{1} = Clim;

end