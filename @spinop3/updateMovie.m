function opts = updateMovie(S, dt, p, opts, t, v, compGrid, plotGrid)
%UPDATEMOVIE   Update the movie when solving a PDE specified by a SPINOP3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
nVars = S.numVars;
dom = S.domain;
Clim = opts{1};
dataToPlot = opts{3};
xx = compGrid{1};
yy = compGrid{2};
zz = compGrid{3};
N = size(xx, 1) - 1;
xxx = plotGrid{1};
yyy = plotGrid{2};
zzz = plotGrid{3};
Nplot = size(xxx, 1) - 1;

% Points for slices:
ttx = trigpts(Nplot, dom(1:2));
tty = trigpts(Nplot, dom(3:4));
ttz = trigpts(Nplot, dom(5:6));
ttx = [ttx; 2*ttx(end) - ttx(end-1)];
tty = [tty; 2*tty(end) - tty(end-1)];
ttz = [ttz; 2*ttz(end) - ttz(end-1)];

% Loop over the variables:
for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vv = dataToPlot(v(idx:idx+N-1,:,:));
    vv = [vv, vv(:,1,:)]; %#ok<*AGROW> add repeated values (periodic endpoints)
    vv = [vv; vv(1,:,:)];
    vv = cat(3, vv, vv(:,:,1));
    
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
        vvv = interp3(xx, yy, zz, vv, xxx, yyy, zzz, 'spline');
    else
        vvv = vv;
    end
    
    % Loop over the surfaces:
    for l = 1:length(p{1,k})
        pkl = p{1,k}(l);
        Sx = ( pkl.XData(1,1) == pkl.XData(1,2) ) && ...
            ( pkl.XData(1,1) == pkl.XData(2,1) ) && ...
            ( pkl.XData(1,1) == pkl.XData(2,2) );
        Sy = ( pkl.YData(1,1) == pkl.YData(1,2) ) && ...
            ( pkl.YData(1,1) == pkl.YData(2,1) ) && ...
            ( pkl.YData(1,1) == pkl.YData(2,2) );
        Sz = ( pkl.ZData(1,1) == pkl.ZData(1,2) ) && ...
            ( pkl.ZData(1,1) == pkl.ZData(2,1) ) && ...
            ( pkl.ZData(1,1) == pkl.ZData(2,2) );
        if ( Sx == 1 )
            pos = p{1,k}(l).XData;
            pos = pos(1);
            [~, id] = min(abs(ttx - pos));
            set(p{1,k}(l), 'cdata', squeeze(vvv(:,id,:)))
        elseif ( Sy == 1 )
            pos = p{1,k}(l).YData;
            pos = pos(1);
            [~, id] = min(abs(tty - pos));
            set(p{1,k}(l), 'cdata', squeeze(vvv(id,:,:)))
        elseif ( Sz == 1 )
            pos = p{1,k}(l).ZData;
            pos = pos(1);
            [~, id] = min(abs(ttz - pos));
            set(p{1,k}(l), 'cdata', squeeze(vvv(:,:,id)))
        end
    end
    ax = p{1,k}.Parent;
    set(ax, 'clim', [Clim(2*(k-1) + 1), Clim(2*(k-1) + 2)])
    
    % Update each title:
    titleString = sprintf('Nx = Ny = Nz = %i (DoFs = %i), dt = %1.1e, t = %.4f', ...
        N, nVars*N^3, dt, t);
    set(p{2,k}, 'String', titleString)
    drawnow
    
end

% Update outputs:
opts{1} = Clim;

end