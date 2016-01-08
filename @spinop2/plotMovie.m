function plotOptions = plotMovie(S, dt, p, plotOptions, t, v, gridPoints)
%PLOTMOVIE   Plot a movie when solving a PDE specified by a SPINOP2.
%   PLOTMOVIE

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
dom = S.domain;
xx = gridPoints{1};
yy = gridPoints{2};
N = size(xx, 1);
nVars = S.numVars;
dataToPlot = plotOptions{2};

% Loop over the variables:
xxplot = [xx, 2*xx(:,end) - xx(:,end-1)];
xxplot = [xxplot; xxplot(1,:)];
yyplot = [yy; 2*yy(end,:) - yy(end-1,:)];
yyplot = [yyplot, yyplot(:,1)];
for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vplot = dataToPlot(v(idx:idx+N-1,:));
    vplot = [vplot, vplot(:,1)]; %#ok<*AGROW>
    vplot = [vplot; vplot(1,:)];
    
    % Update it:
    set(p{k}, 'xdata', xxplot, 'ydata', yyplot, 'zdata', vplot)
    set(p{k}.Parent, 'xlim', [dom(1), dom(2)], 'ylim', [dom(3) dom(4)])

    % Update title:
    if ( k == 1 )
        title(p{k}.Parent, sprintf(['N = %i (DoFs = %i), dt = %1.1e, ', ...
            't = %.4f'], N, nVars*N^2, dt, t))
    end
    drawnow
end

end