function [p, opts] = initializeMovie(S, dt, pref, v, compGrid, plotGrid)
%INITIALIZEMOVIE   Initialize a movie when solving a PDE specified by a SPINOP.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Note: P is a NVARSx2 CELL-ARRAY that stores the NVARS plots in the first 
% row, and the NVARS titles in the second row. OPTS is a 2x1 CELL-ARRAY that 
% stores the limits of the y-axis in OPTS{1} and what kind of data to plot in 
% OPTS{2} (real/imag/abs; when the data is complex-valued).

% Set-up:
dom = S.domain;                       % Spatial domain
nVars = S.numVars;                    % Number of variables (>1 for systems)
vscale = max(abs(v));                 % Scale of the solution 
dataplot = str2func(pref.dataplot);   % Plot 'abs', 'real' or 'imag'
xx = compGrid{1};                     % Computation grid
xxx = plotGrid{1};                    % Movie grid
N = size(xx, 1) - 1;                  % Size of computation grid
Nplot = size(xxx, 1) - 1;             % Size of movie grid 
FS = 'fontsize'; fs = 12;             % Fontsize for title 

% Loop over the variables:
p = cell(2, nVars); clf reset
for k = 1:nVars
    
    % Extract each variable:
    idx = (k-1)*N + 1;
    vv = dataplot(v(idx:idx+N-1));
    vv = [vv; vv(1)]; %#ok<*AGROW> add repeated values (periodic endpoints)
    
    % Get the YLIM for the y-axis:
    if ( isempty(pref.Ylim) == 1 )
        Ylim(2*(k-1) + 1) = min(vv) - .1*vscale;
        Ylim(2*(k-1) + 2) = max(vv) + .1*vscale;
    else
        Ylim = pref.Ylim;
    end
    
    % Interpolate each variable on a finer grid:
    if ( Nplot > N )
        vvv = interp1(xx, vv, xxx, 'spline');
    else
        vvv = vv;
    end
    
    % Plot each variable:
    subplot(1, nVars, k)
    p{1,k} = plot(xxx, vvv, 'linewidth', 3);
    axis([dom(1), dom(2), Ylim(2*(k-1) + 1), Ylim(2*(k-1) + 2)])
    if ( nVars == 1 )
        xlabel('x'), ylabel('u(t,x)'), grid on
    else
        xlabel('x'), ylabel(['u_',num2str(k),'(t,x)']), grid on
    end
    set(gca, FS, fs), box on
    
    % Plot each title:
    titleString = sprintf('N = %i (DoFs = %i), dt = %1.1e, t = %.4f', N, ...
        nVars*N, dt, 0);
    p{2,k} = title(titleString);
    drawnow
    
end

% Ask the user to press SPACE:
state = pause;
if ( strcmpi(state, 'on') == 1 )
    disp('Type <space> when ready.')
end
shg, pause

% Outputs:
opts{1} = Ylim;
opts{2} = dataplot;

end