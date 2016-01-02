function [p, plotOption] = initializeMovie(S, dt, pref, v, gridPoints)
%INITIALIZEMOVIE   Initialize a movie when solving a PDE specified by a SPINOP2.
%   INITIALIZEMOVIE

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
dom = S.domain;
nVars = S.numVars;
viewSpec = pref.view;
defaultPref = spinpref2();
defaultView = defaultPref.view;
while ( length(viewSpec) < 2*nVars )
    viewSpec = [viewSpec, defaultView];
end
xx = gridPoints{1};
yy = gridPoints{2};
N = size(xx, 1);

% Plot:
xxplot = [xx, 2*xx(:,end) - xx(:,end-1)];
xxplot = [xxplot; xxplot(1,:)];
yyplot = [yy; 2*yy(end,:) - yy(end-1,:)];
yyplot = [yyplot, yyplot(:,1)];
p = cell(nVars, 1); clf
for k = 1:nVars
    idx = (k-1)*N + 1;
    vplot = v(idx:idx+N-1,:);
    vplot = [vplot, vplot(:,1)]; %#ok<*AGROW>
    vplot = [vplot; vplot(1,:)];
    subplot(nVars, 1, k)
    p{k} = surf(xxplot, yyplot, vplot, 'edgecolor', 'none');
    axis([dom(1) dom(2) dom(3) dom(4)])
    view(viewSpec(2*(k - 1) + 1 : 2*(k - 1) + 2)), colorbar
    xlabel('x'), ylabel('y'), set(gca, 'FontSize', 16), box on
    if ( k == 1 )
        title(p{k}.Parent, sprintf('N = %i, dt = %1.1e, t = %.4f', N, dt, 0))
    end
    drawnow
end
disp('Type <enter> when ready.'), shg, pause
plotOption = viewSpec;

end