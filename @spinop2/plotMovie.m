function options = plotMovie(S, dt, p, options, t, v, dataGrid, plotGrid)
%PLOTMOVIE   Plot a movie when solving a PDE specified by a SPINOP2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
dom = S.domain;
nVars = S.numVars;
Clim = options{1};
dataToPlot = options{3};
xx = dataGrid{1};
yy = dataGrid{2};
N = size(xx, 1) - 1;
xxx = plotGrid{1};
yyy = plotGrid{2};

for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vv = dataToPlot(v(idx:idx+N-1,:));
    vv = [vv, vv(:,1)]; %#ok<*AGROW>
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
    vvv = interp2(xx, yy, vv, xxx, yyy, 'spline');
    
    % Update each variable:
    set(p{k}, 'xdata', xxx, 'ydata', yyy, 'zdata', vvv)
    set(p{k}.Parent, 'xlim', [dom(1), dom(2)], 'ylim', [dom(3) dom(4)])
    set(p{k}.Parent, 'clim', [Clim(2*(k-1) + 1), Clim(2*(k-1) + 2)])
    drawnow
    
end

% Update title:
titleString = sprintf('Nx = Ny = %i (DoFs = %i), dt = %1.1e, t = %.4f', N, ...
    nVars*N^2, dt, t);
set(p{nVars + 1}, 'String', titleString)

% Update outputs:
options{1} = Clim;

end