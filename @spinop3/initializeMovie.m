function [p, plotOption] = initializeMovie(S, dt, pref, v, gridPoints)
%INITIALIZEMOVIE   Initialize a movie when solving a PDE specified by a SPINOP3.
%   INITIALIZEMOVIE

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
xx = gridPoints{1};
yy = gridPoints{2};
zz = gridPoints{3};
N = size(xx, 1);
nVars = S.numVars;
dom = S.domain;
xxplot = [xx, 2*xx(:,end,:) - xx(:,end-1,:)];
xxplot =  [xxplot; xxplot(1,:,:)];
xxplot = cat(3, xxplot, xxplot(:,:,1));
yyplot = [yy; 2*yy(end,:,:) - yy(end-1,:,:)];
yyplot = [yyplot, yyplot(:,1,:)];
yyplot = cat(3, yyplot, yyplot(:,:,1));
zzplot = cat(3, zz, 2*zz(:,:,end) - zz(:,:,end-1));
zzplot = [zzplot; zzplot(1,:,:)];
zzplot = [zzplot, zzplot(:,1,:)];

% Loop over the variables:
clf
for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vplot = real(v(idx:idx+N-1,:,:));
    vplot = [vplot, vplot(:,1,:)]; %#ok<*AGROW>
    vplot = [vplot; vplot(1,:,:)];
    vplot = cat(3, vplot, vplot(:,:,1));
    
    % Plot it:
    subplot(nVars, 1, k)
    isosurface(xxplot, yyplot, zzplot, vplot)
    axis([dom(1) dom(2) dom(3) dom(4) dom(5) dom(6)]), colorbar, camlight 
    xlabel('x'), ylabel('y'), zlabel('z'), set(gca, 'FontSize', 16), box on
    if ( k == 1 )
        title(sprintf('N = %i, dt = %1.1e, t = %.4f', N, dt, 0))
    end
    drawnow
    
end
disp('Type <space> when ready.'), pause
p = [];
plotOption = [];

end