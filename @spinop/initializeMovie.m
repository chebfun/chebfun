function [p, plotOption] = initializeMovie(S, dt, pref, v, gridPoints)
%INITIALIZEMOVIE   Initialize a movie when solving a PDE specified by a SPINOP.
%   INITIALIZEMOVIE

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
dom = S.domain;
xx = gridPoints;
N = size(xx, 1);
nVars = S.numVars;

% Plot:
xxplot = [xx; 2*xx(end) - xx(end-1)];
vscale = max(abs(v));
p = cell(nVars, 1); clf
for k = 1:nVars
    idx = (k-1)*N + 1;
    vplot = v(idx:idx+N-1);
    vplot = [vplot; vplot(1)]; %#ok<*AGROW>
    if ( isempty(pref.Ylim) == 1 )
        Ylim(2*(k-1) + 1) = min(vplot) - .1*vscale;
        Ylim(2*(k-1) + 2) = max(vplot) + .1*vscale;
    else
        Ylim = pref.Ylim;
    end
    subplot(nVars, 1, k)
    p{k} = plot(xxplot, vplot, 'linewidth', 3);
    axis([dom(1), dom(2), Ylim(2*(k-1) + 1), Ylim(2*(k-1) + 2)])
    if ( nVars == 1 )
        xlabel('x'), ylabel('u(t,x)'), grid on
    else
        xlabel('x'), ylabel(['u_',num2str(k),'(t,x)']), grid on
    end
    set(gca, 'FontSize', 16), box on
    if ( k == 1 )
        title(p{k}.Parent, sprintf('N = %i, dt = %1.1e, t = %.4f', N, dt, 0))
    end
    drawnow
end
disp('Type <enter> when ready.'), shg, pause
plotOption = Ylim;

end