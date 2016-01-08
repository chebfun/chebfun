function [p, plotOptions] = initializeMovie(S, dt, pref, v, gridPoints)
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
dataToPlot = str2func(pref.dataToPlot);
dom = S.domain;
tt = trigpts(N, dom);
tt = [tt; 2*tt(end) - tt(end-1)];
if ( isempty(pref.slices) == 1 )
    leftPos = floor(N/10);
    midPos = floor(N/2) + 1;
    rightPos = floor(9*N/10);
    Sx = tt(midPos);
    Sy = Sx;
    Sz = [tt(leftPos); tt(rightPos)];
else
    slices = pref.slices;
    Sx = slices{1};
    for k = 1:length(Sx)
        pos = Sx(k);
        [~, id] = min(abs(tt-pos));
        Sx(k) = tt(id);
    end
    Sy = slices{2};
    for k = 1:length(Sy)
        pos = Sy(k);
        [~, id] = min(abs(tt-pos));
        Sy(k) = tt(id);
    end
    Sz = slices{3};
    for k = 1:length(Sz)
        pos = Sz(k);
        [~, id] = min(abs(tt-pos));
        Sz(k) = tt(id);
    end
end
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
p = cell(nVars, 1); clf
for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vplot = dataToPlot(v(idx:idx+N-1,:,:));
    vplot = [vplot, vplot(:,1,:)]; %#ok<*AGROW>
    vplot = [vplot; vplot(1,:,:)];
    vplot = cat(3, vplot, vplot(:,:,1));
    
    % Plot it:
    subplot(nVars, 1, k) 
    p{k} = slice(xxplot, yyplot, zzplot, vplot, Sx, Sy, Sz);
    axis([dom(1) dom(2) dom(3) dom(4) dom(5) dom(6)]), colorbar
    xlabel('x'), ylabel('y'), zlabel('z'), set(gca, 'FontSize', 16), box on
    if ( k == 1 )
        title(sprintf(['N = %i (DoFs = %i), dt = %1.1e, ', 't = %.4f'], N, ...
            nVars*N^3, dt, 0))
    end
    drawnow
    
end
disp('Type <space> when ready.'), shg, pause
plotOptions{1} = {Sx, Sy, Sz};
plotOptions{2} = dataToPlot;

end