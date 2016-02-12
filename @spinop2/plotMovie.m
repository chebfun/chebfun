function plotOptions = plotMovie(S, dt, p, plotOptions, t, v, gridPoints)
%PLOTMOVIE   Plot a movie when solving a PDE specified by a SPINOP2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
dom = S.domain;
xx = gridPoints{1};
yy = gridPoints{2};
N = size(xx, 1);
nVars = S.numVars;
Clim = plotOptions{1};
dataToPlot = plotOptions{3};

% Grid of the computation:
xxplot = [xx, 2*xx(:,end) - xx(:,end-1)];
xxplot = [xxplot; xxplot(1,:)];
yyplot = [yy; 2*yy(end,:) - yy(end-1,:)];
yyplot = [yyplot, yyplot(:,1)];

% Finer grid for interploation:
Nplot = max(N, 256);
[xxxplot, yyyplot] = meshgrid(trigpts(Nplot, dom));
xxxplot = [xxxplot, 2*xxxplot(:,end) - xxxplot(:,end-1)];
xxxplot = [xxxplot; xxxplot(1,:)];
yyyplot = [yyyplot; 2*yyyplot(end,:) - yyyplot(end-1,:)];
yyyplot = [yyyplot, yyyplot(:,1)];

for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vvplot = dataToPlot(v(idx:idx+N-1,:));
    vvplot = [vvplot, vvplot(:,1)]; %#ok<*AGROW>
    vvplot = [vvplot; vvplot(1,:)];
    
    % Change axes if necessary:
    if ( nargout == 1 )
        minvnew = min(vvplot(:));
        maxvnew = max(vvplot(:));
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
    vvvplot = interp2(xxplot, yyplot, vvplot, xxxplot, yyyplot, 'spline');
    
    % Update each variable:
    set(p{k}, 'xdata', xxxplot, 'ydata', yyyplot, 'zdata', vvvplot)
    set(p{k}.Parent, 'xlim', [dom(1), dom(2)], 'ylim', [dom(3) dom(4)])
    set(p{k}.Parent, 'clim', [Clim(2*(k-1) + 1), Clim(2*(k-1) + 2)])
    drawnow
    
end

% Update title:
titleString = sprintf('Nx = Ny = %i (DoFs = %i), dt = %1.1e, t = %.4f', N, ...
    nVars*N^2, dt, t);
set(p{nVars + 1}, 'String', titleString)

% Update outputs:
plotOptions{1} = Clim;

end