function [p, plotOptions] = initializeMovie(S, dt, pref, v, gridPoints)
%INITIALIZEMOVIE   Initialize a movie when solving a PDE specified by a SPINOP2.
%   INITIALIZEMOVIE

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
dom = S.domain;
nVars = S.numVars;
viewSpec = pref.view;
dataToPlot = str2func(pref.dataToPlot);
defaultPref = spinpref2();
defaultView = defaultPref.view;
while ( length(viewSpec) < 2*nVars )
    viewSpec = [viewSpec, defaultView];
end
xx = gridPoints{1};
yy = gridPoints{2};
N = size(xx, 1);
xxplot = [xx, 2*xx(:,end) - xx(:,end-1)];
xxplot = [xxplot; xxplot(1,:)];
yyplot = [yy; 2*yy(end,:) - yy(end-1,:)];
yyplot = [yyplot, yyplot(:,1)];

% Loop over the variables:
p = cell(nVars, 1); clf
for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vplot = dataToPlot(v(idx:idx+N-1,:));
    vplot = [vplot, vplot(:,1)]; %#ok<*AGROW>
    vplot = [vplot; vplot(1,:)];
    
    % Plot it:
    subplot(nVars, 1, k)
    p{k} = surf(xxplot, yyplot, vplot, 'edgecolor', 'none');
    axis([dom(1) dom(2) dom(3) dom(4)])
    view(viewSpec(2*(k - 1) + 1 : 2*(k - 1) + 2)), colorbar
    xlabel('x'), ylabel('y'), set(gca, 'FontSize', 16), box on
    
    % Title:
    if ( k == 1 )
        title(p{k}.Parent, sprintf(['N = %i (DoFs = %i), dt = %1.1e, ', ...
            't = %.4f'], N, nVars*N^2, dt, 0))
    end
    drawnow
    
end
disp('Type <space> when ready.'), shg, pause
plotOptions{1} = viewSpec;
plotOptions{2} = dataToPlot;

end