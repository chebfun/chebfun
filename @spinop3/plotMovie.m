function plotOptions = plotMovie(S, dt, p, plotOptions, t, v, gridPoints)
%PLOTMOVIE   Plot a movie when solving a PDE specified by a SPINOP3.
%   PLOTMOVIE

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
xx = gridPoints{1};
yy = gridPoints{2};
zz = gridPoints{3};
N = size(xx, 1);
nVars = S.numVars;
dataToPlot = plotOptions{2};
dom = S.domain;

% Grid of the computation:
xxplot = [xx, 2*xx(:,end,:) - xx(:,end-1,:)];
xxplot =  [xxplot; xxplot(1,:,:)];
xxplot = cat(3, xxplot, xxplot(:,:,1));
yyplot = [yy; 2*yy(end,:,:) - yy(end-1,:,:)];
yyplot = [yyplot, yyplot(:,1,:)];
yyplot = cat(3, yyplot, yyplot(:,:,1));
zzplot = cat(3, zz, 2*zz(:,:,end) - zz(:,:,end-1));
zzplot = [zzplot; zzplot(1,:,:)];
zzplot = [zzplot, zzplot(:,1,:)];

% Finer grid for interploation:
Nplot = max(N, 64);
tt = trigpts(Nplot, dom);
[xxxplot, yyyplot, zzzplot] = meshgrid(trigpts(Nplot, dom));
xxxplot = [xxxplot, 2*xxxplot(:,end,:) - xxxplot(:,end-1,:)];
xxxplot =  [xxxplot; xxxplot(1,:,:)];
xxxplot = cat(3, xxxplot, xxxplot(:,:,1));
yyyplot = [yyyplot; 2*yyyplot(end,:,:) - yyyplot(end-1,:,:)];
yyyplot = [yyyplot, yyyplot(:,1,:)];
yyyplot = cat(3, yyyplot, yyyplot(:,:,1));
zzzplot = cat(3, zzzplot, 2*zzzplot(:,:,end) - zzzplot(:,:,end-1));
zzzplot = [zzzplot; zzzplot(1,:,:)];
zzzplot = [zzzplot, zzzplot(:,1,:)];

% Loop over the variables:
for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vvplot = dataToPlot(v(idx:idx+N-1,:,:));
    vvplot = [vvplot, vvplot(:,1,:)]; %#ok<*AGROW>
    vvplot = [vvplot; vvplot(1,:,:)];
    vvplot = cat(3, vvplot, vvplot(:,:,1));
    
    % Interpolate each variable on a finer grid:
    vvvplot = interp3(xxplot, yyplot, zzplot, vvplot, xxxplot, yyyplot, ...
        zzzplot, 'spline');
    
    % Loop over the surfaces:
    subplot(nVars, 1, k)
    for l = 1:length(p{k})
        pkl = p{k}(l);
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
            pos = p{k}(l).XData;
            pos = pos(1);
            [~, id] = min(abs(tt-pos));
            set(p{k}(l), 'xdata', squeeze(xxxplot(:,id,:)))
            set(p{k}(l), 'ydata', squeeze(yyyplot(:,id,:)))
            set(p{k}(l), 'zdata', squeeze(zzzplot(:,id,:)))
            set(p{k}(l), 'cdata', squeeze(vvvplot(:,id,:)))
        elseif ( Sy == 1 )
            pos = p{k}(l).YData;
            pos = pos(1);
            [~, id] = min(abs(tt-pos));
            set(p{k}(l), 'xdata', squeeze(xxxplot(id,:,:)))
            set(p{k}(l), 'ydata', squeeze(yyyplot(id,:,:)))
            set(p{k}(l), 'zdata', squeeze(zzzplot(id,:,:)))
            set(p{k}(l), 'cdata', squeeze(vvvplot(id,:,:)))
        elseif ( Sz == 1 )
            pos = p{k}(l).ZData;
            pos = pos(1);
            [~, id] = min(abs(tt-pos));
            set(p{k}(l), 'xdata', squeeze(xxxplot(:,:,id)))
            set(p{k}(l), 'ydata', squeeze(yyyplot(:,:,id)))
            set(p{k}(l), 'zdata', squeeze(zzzplot(:,:,id)))
            set(p{k}(l), 'cdata', squeeze(vvvplot(:,:,id)))
        end
    end
    drawnow
    
end

% Update title:
titleString = sprintf('Nx = Ny = Nz = %i (DoFs = %i), dt = %1.1e, t = %.4f', ...
    N, nVars*N^3, dt, t);
set(p{nVars + 1}, 'String', titleString)

end