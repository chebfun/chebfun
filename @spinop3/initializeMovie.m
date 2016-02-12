function [p, plotOptions] = initializeMovie(S, dt, pref, v, gridPoints)
%INITIALIZEMOVIE   Initialize a movie when solving a PDE specified by a SPINOP3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
xx = gridPoints{1};
yy = gridPoints{2};
zz = gridPoints{3};
N = size(xx, 1);
nVars = S.numVars;
vscale = max(abs(v(:)));
dataToPlot = str2func(pref.dataToPlot);
dom = S.domain;
tt = trigpts(N, dom);
tt = [tt; 2*tt(end) - tt(end-1)];
if ( isempty(pref.slices) == 1 )
    Sx = tt(floor(N/2) + 1);
    Sy = Sx;
    Sz = Sx;
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
Nplot = max(N, 100);
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
p = cell(nVars + 1, 1); clf reset
for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vvplot = dataToPlot(v(idx:idx+N-1,:,:));
    vvplot = [vvplot, vvplot(:,1,:)]; %#ok<*AGROW>
    vvplot = [vvplot; vvplot(1,:,:)];
    vvplot = cat(3, vvplot, vvplot(:,:,1));
    
    % Get the CLIM for the colorbar:
    if ( isempty(pref.Clim) == 1 )
        Clim(2*(k-1) + 1) = min(vvplot(:)) - .1*vscale;
        Clim(2*(k-1) + 2) = max(vvplot(:)) + .1*vscale;
    else
        Clim = pref.Clim;
    end
    
    % Interpolate each variable on a finer grid:
    vvvplot = interp3(xxplot, yyplot, zzplot, vvplot, xxxplot, yyyplot, ...
        zzzplot, 'spline');
    
    % Plot each variable:
    subplot(1, nVars, k) 
    p{k} = slice(xxxplot, yyyplot, zzzplot, vvvplot, Sx, Sy, Sz);
    set(p{k}, 'edgecolor', 'none')
    ax = p{k}.Parent; set(ax, 'clim', [Clim(2*(k-1) + 1), Clim(2*(k-1) + 2)])
    axis([dom(1) dom(2) dom(3) dom(4) dom(5) dom(6)]), colorbar
    xlabel('x'), ylabel('y'), zlabel('z'), set(gca, 'FontSize', 16), box on
    drawnow
    
end

% Title:
titleString = sprintf('Nx = Ny = Nz = %i (DoFs = %i), dt = %1.1e, t = %.4f', ...
    N, nVars*N^3, dt, 0);
set(gcf, 'NextPlot', 'add');
ax = axes;
h = title(titleString);
set(ax, 'Visible', 'off', 'HandleVisibility', 'off', 'Fontsize', 16)
set(h, 'Visible', 'on', 'Position', [.47 1.01 .5])

% Ask the user to press SPACE:
state = pause;
if ( strcmpi(state, 'on') == 1 )
    disp('Type <space> when ready.')
end
shg, pause

% Outputs:
p{nVars + 1} = h;
plotOptions{1} = Clim;
plotOptions{2} = {Sx, Sy, Sz};
plotOptions{3} = dataToPlot;

end