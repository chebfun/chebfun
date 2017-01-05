function [uOut, tOut, computingTime] = solvepde(varargin)
%SOLVEPDE   Solve a PDE defined by a SPINOP, SPINOP2, SPINOP3 or SPINOPSPHERE.
%   SOLVEPDE is called by SPIN, SPIN2, SPIN3 and SPINSPHERE. It is not called
%   directly by the user. Appropriate help texts can be found in SPIN, SPIN2,
%   SPIN3 and SPINSPHERE.
%
% See also SPIN, SPIN2, SPIN3, SPINSPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% Parse inputs:

% SOLVEPDE has been called by SPIN/SPIN2/SPIN3. The inputs have been parsed in
% those files and are expeceted to be:
%
% OPTION 1.     SOLVEPDE(S, N, DT), S is a SPINOPERATOR object, N is the number
%               of grid points and DT is the time-step.
%
% OPTION 2.     SOLVEPDE(S, N, DT, PREF), PREF is a SPINPREFERENCE.

% Get the inputs:
pref = [];
if ( nargin == 3 ) % OPTION 1
    S = varargin{1};
    N = varargin{2};
    dt = varargin{3};
elseif ( nargin == 4 ) % OPTION 2
    S = varargin{1};
    N = varargin{2};
    dt = varargin{3};
    pref = varargin{4};
end

% Dimension:
dim = getDimension(S);

% Create a SPINPREFERENCE object if none:
if ( isempty(pref) == 1 )
    pref = getPreference(S);
end

%% Pre-processing:

% Time interval TSPAN:
tspan = S.tspan;

% Initial condition U0:
u0 = S.init;
if ( isa(u0, 'chebfun') == 1 || isa(u0, 'chebfun2') == 1 || ...
        isa(u0, 'chebfun3') == 1 || isa(u0, 'spherefun') == 1 )
    u0 = chebmatrix(u0);
elseif ( isa(u0, 'chebfun2v') == 1 || isa(u0, 'chebfun3v') == 1 || ...
        isa(u0, 'spherefunv') == 1 )
    temp = chebmatrix(u0(1));
    for k = 2:size(u0, 1)
        temp(k,1) = u0(k);
    end
    u0 = temp;
end

% Number of variables:
nVars = S.numVars;

% Space interval DOM and final time TF:
dom = u0{1}.domain;
tf = tspan(end);

% Throw an error if the domains of the initial condition and the SPINOPERATOR
% are different:
if ( isempty(S) == 0 )
    if ( isequal(S.domain, dom) == 0 )
        error('SPINOPERATOR:solvepde', ['The initial condition and the ', ...
            'operator do not live on the same domain.'])
    end
end

% Get the preferences:
dealias = pref.dealias;   % Use a dealiasing procedure if DEALIAS = 1
iterplot = pref.iterplot; % plot every ITERPLOT iterations if 'movie'
plotStyle = pref.plot;    % Plotting options

% Create a time-stepping scheme:
schemeName = pref.scheme;
K = spinscheme(schemeName);
q = K.steps;  % Number of steps of the scheme (q>1 for multistep schemes)

% Diagonal SPINOPERATOR objects (1D/2D/3D) use exponential integtrators while
% nondiagonal SPINOPERATOR objects (sphere) use IMEX schemes. Check we're using
% the right scheme:
if ( isDiag(S) == 1 ) % 1D/2D/3D
    if ( strcmpi(K.type, 'expint') ~= 1 )
        error(['Use exponential integrators with SPIN/SPIN2/SPIN3. ', ...
            'See HELP/SPINSCHEME.'])
    end
else % sphere
    if ( strcmpi(K.type, 'imex') ~= 1 )
        error('Use IMEX schemes with SPINSPHERE. See HELP/SPINSCHEME.')
    end
end

% Exponential integrators use contour integrals for computing the phi-functions:
if ( strcmpi(K.type, 'expint') == 1 )
    M = pref.M;                       % Number of points for complex means
end

% Operators: linear part L, and nonlinear parts Nc (in coefficient space) and 
% Nv (in value space). The nonlinear opeartor, acting on the Fourier coeffs
% is written as 
%           
%        N(coeffs) = Nc(fft(Nv(ifft(coeffs)))).
%          
% For example, for u_t = Lu + (u^2)_x,
%
%        N(coeffs) = 1i*fft(ifft(coeffs)^2) with Nv(u) = u^2 and Nc(u) = 1i*u.
%
[L, Nc] = discretize(S, N);
Nv = S.nonlinearPartVals;

% Set-up spatial grid:
grid = getGrid(S, N, dom);

% Get the values to coeffs and coeffs to values transform:
v2c = getVals2CoeffsTransform(S);
c2v = getCoeffs2ValsTransform(S);

% Get the values VINIT and Fourier coeffs CINIT of the initial condition.
% Also store the nonlinear evaluation of the initial Fourier coefficients 
% in NCINIT.
if ( isDiag(S) == 1 ) % the linear part of the operartor is diagonal (1D/2D/3D)
    
    % Initial values VINIT:
    vInit = [];
    for k = 1:nVars
        if ( dim == 1 ) % 1D
            xx = grid{1};
            vInit = [vInit; feval(u0{k}, xx)]; %#ok<*AGROW>
        elseif ( dim == 2 ) % 2D 
            xx = grid{1};
            yy = grid{2};
            vInit = [vInit; feval(u0{k}, xx, yy)];
        elseif ( dim == 3 ) % 3D
            xx = grid{1}; 
            yy = grid{2};
            zz = grid{3};
            vInit = [vInit; feval(u0{k}, xx, yy, zz)];
        end
    end
    
    % Initial Fourier coefficients CINIT:
    cInit{1} = [];
    for k = 1:nVars
        idx = (k-1)*N + 1;
        cInit{1} = [cInit{1}; v2c(vInit(idx:idx+N-1,:,:))]; 
    end
    
    % Nonlinear evaluation of the initial Fourier coefficients NCINIT:
    vals = Nv(vInit);            % apply NV in value space
    coeffs = v2c(vals(1:N,:,:)); % transform to coefficients
    for k = 1:nVars-1
        idx = k*N + 1;
        coeffs = [coeffs; v2c(vals(idx:idx+N-1,:,:))];
    end
    coeffs = Nc.*coeffs;         % apply NC in coefficient space
    NcInit{1} = coeffs;
    
else % the linear part of the operartor is not diagonal (sphere)
    
    if ( strcmpi(dim, 'unit sphere') == 1 ) % sphere
        
        % Initial Fourier coefficients CINIT:
        cInit{1} = reshape(coeffs2(u0{1}, N, N), N*N, 1);
        for k = 2:nVars
            cInit{1} = [cInit{1}; reshape(coeffs2(u0{k}, N, N), N*N, 1);];
        end
        
        % Initial values VINIT:
        vInit = c2v(cInit{1}(1:N^2));
        for k = 1:nVars-1
            idx = k*N^2 + 1;
            vInit = [vInit; c2v(cInit{1}(idx:idx+N^2-1))];
        end
        
        % Nonlinear evaluation of the initial Fourier coefficients NCINIT:
        vals = Nv(vInit);           % apply NV in value space
        coeffs = v2c(vals(1:N,:));  % transform to coefficients
        for k = 1:nVars-1
            idx = k*N + 1;
            coeffs = [coeffs; v2c(vals(idx:idx+N-1,:))];
        end
        coeffs = Nc*coeffs;         % apply NC in coefficient space
        NcInit{1} = coeffs;
    end
end

% Get enough initial data when using a multistep scheme:
if ( q > 1 )
    if ( strcmpi(K.type, 'expint') == 1 ) % exponential integrators
        [cInit, NcInit] = startMultistep(K, dt, L, Nc, Nv, pref, S, cInit, ...
            NcInit);
    else % imex
        % [TODO]: Implement multistep IMEX schemes.
        error('Multistep IMEX schemes are not supported.')
    end
end

% Compute the coefficients of the scheme for the exponential integrators:
if ( strcmpi(K.type, 'expint') == 1 )
    schemeCoeffs = computeCoeffs(K, dt, L, M, S);
else
    schemeCoeffs = [];
end

% Indexes for dealiasing:
toOne = floor(N/2) + 1 - ceil(N/6):floor(N/2) + ceil(N/6);
if ( dim == 1 )
    ind = false(N, 1);
    ind(toOne) = 1;
elseif ( dim == 2 )
    ind = false(N, N);
    ind(toOne, toOne) = 1;
elseif ( dim == 3 )
    ind = false(N, N, N);
    ind(toOne, toOne, toOne) = 1;
end
ind = repmat(ind, nVars, 1);

% Values VOUT and times TOUT to output:
vOut{1} = vInit;
tOut(1) = 0;

% Values VWATER to plot if using WATERFALL:
if ( strcmpi(plotStyle, 'waterfall') == 1 )
    vWater{1} = vInit;
    twater = 0;
end

% Create grids for plotting, and plot initial condition if using MOVIE:
if ( strcmpi(plotStyle, 'movie') == 1 )
    Nplot = max(N, pref.Nplot);
    if ( dim == 1 )
        computationGrid = {xx};
        plottingGrid = {trigpts(Nplot, dom)};
    elseif ( dim == 2 || strcmpi(dim, 'unit sphere') == 1 )
        computationGrid = {xx; yy};
        ttx = trigpts(Nplot, dom(1:2));
        tty = trigpts(Nplot, dom(3:4));
        [xxx, yyy] = meshgrid(ttx, tty);
        plottingGrid = {xxx; yyy};
    elseif ( dim == 3 );
        computationGrid = {xx; yy; zz};
        ttx = trigpts(Nplot, dom(1:2));
        tty = trigpts(Nplot, dom(3:4));
        ttz = trigpts(Nplot, dom(5:6));
        [xxx, yyy, zzz] = meshgrid(ttx, tty, ttz);
        plottingGrid = {xxx; yyy; zzz};
    end
    computationGrid = reshapeGrid(S, computationGrid); 
    plottingGrid = reshapeGrid(S, plottingGrid); 
    [p, options] = initializeMovie(S, dt, pref, vInit, computationGrid, plottingGrid);
end

%% Time-stepping loop:

tic
iter = 0;
valuesUpdated = 0;
t = (q-1)*dt;
pos = 2;
cOld = cInit;
NcOld = NcInit;
while ( t < tf )
    
    % One step in time with DT and N points:
    [cNew, NcNew] = oneStep(K, schemeCoeffs, Nc, Nv, nVars, cOld, NcOld);
    if ( any(isnan(cNew{1})) )
        error('The solution blew up. Try a smaller time-step.')
    end
    iter = iter + 1;
    t = (iter + q - 1)*dt;
    cOld = cNew;
    NcOld = NcNew;
    
    % Dealiasing procedure:
    if ( strcmpi(dealias, 'on') == 1 )
        cNew{1}(ind) = 0;
    end
    
    % Plot every ITERPLOT iterations if using MOVIE:
    if ( strcmpi(plotStyle, 'movie') == 1 && mod(iter,iterplot) == 0 )
        v = [];
        for k = 1:nVars
            idx = (k-1)*N + 1;
            temp = ifftn(cNew{1}(idx:idx+N-1,:,:));
            if ( max(abs(imag(temp(:)))) < max(abs(temp(:)))*1e-10 )
                temp = real(temp);
            end
            v = [v; temp];
        end
        valuesUpdated = 1;
        if ( dim == 1 )
            isLimGiven = ~isempty(pref.Ylim);
        elseif ( dim == 2 || dim == 3 )
            isLimGiven = ~isempty(pref.Clim);
        end
        if ( isLimGiven == 1 )
            plotMovie(S, dt, p, options, t, v, computationGrid, plottingGrid);
        else
            options = plotMovie(S, dt, p, options, t, v, computationGrid, plottingGrid);
        end
        
    % Store the values every ITERPLOT iterations if using WATERFALL:
    % (Remark: Only in dimension 1.)
    elseif ( strcmpi(plotStyle, 'waterfall') == 1 && mod(iter, iterplot) == 0 )
        v = [];
        for k = 1:nVars
            idx = (k-1)*N + 1;
            temp = ifft(cNew{1}(idx:idx+N-1));
            if ( max(abs(imag(temp(:)))) < max(abs(temp(:)))*1e-10 )
                temp = real(temp);
            end
            v = [v; temp];
        end
        valuesUpdated = 1;
        vWater{iter/iterplot + 1} = v;
        twater = [twater, t];
    end
    
    % Output the solution if T correponds to an entry of TSPAN:
    if ( abs(t - tspan(pos)) < 1e-10 )
        if ( valuesUpdated == 0 )
            v = [];
            for k = 1:nVars
                idx = (k-1)*N + 1;
                temp = ifftn(cNew{1}(idx:idx+N-1,:,:));
                if ( max(abs(imag(temp(:)))) < max(abs(temp(:)))*1e-10 )
                    temp = real(temp);
                end
                v = [v; temp];
            end
        end
        vOut{pos} = v;
        tOut(pos) = t;
        pos = pos + 1;
        if ( pos > length(tspan) )
            break
        end
    end
    
end
computingTime = toc;

%% Post-processing:

% Make sure that the solution at TF has been plotted if using MOVIE:
if ( strcmpi(plotStyle, 'movie') == 1 )
    plotMovie(S, dt, p, options, t, v, computationGrid, plottingGrid);
    set(gcf, 'NextPlot', 'replace')
end

% Use WATERFALL if using WATERFALL:
if ( strcmpi(plotStyle, 'waterfall') == 1 )
    clf reset
    for k = 1:nVars
        uwater = [];
        for l = 1:size(vWater, 2)
            N = length(vWater{l})/nVars;
            idx = (k-1)*N + 1;
            uwater = [ uwater, chebfun(real(vWater{l}(idx:idx+N-1)), dom, ...
                'trig') ];
        end
        subplot(1, nVars, k)
        waterfall(uwater, twater), axis([dom(1), dom(end), 0, tf])
        set(gca, 'FontSize', 16), box on
        xlabel('x'), ylabel('t')
        if ( nVars == 1 )
            zlabel('u(t,x)')
        else
            zlabel(['u_',num2str(k),'(t,x)'])
        end
        view([10 70])
    end
end

% Get the right type of CHEBFUN:
if ( dim == 1 )
    fun = @chebfun;
elseif ( dim == 2 )
    fun = @chebfun2;
elseif ( dim == 3 )
    fun = @chebfun3;
    
    % The data come from MESHGRID, need to permute them because the CHEBFUN3
    % constructor assumes that the data come from NDGRID:
    for l = 1:size(vOut, 2)
        for k = 1:nVars
            idx = (k-1)*N + 1;
            vals = vOut{l}(idx:idx+N-1,:,:);
            vOut{l}(idx:idx+N-1,:,:) = permute(vals, [2 1 3]);
        end
    end
    
end

% Output a CHEBFUN/CHEBMATRIX from values VOUT:
if ( length(tspan) == 2 ) % e.g., tspan = [0 T], output only the solution at T
    N = length(vOut{end})/nVars;
    if ( nVars == 1 )
        uOut = fun(vOut{end}(1:N,:,:), dom, 'trig');
    else % for systems, use a CHEBAMTRIX
        uOut = chebmatrix(fun(vOut{end}(1:N,:,:), dom, 'trig'));
        for k = 2:nVars
            idx = (k-1)*N + 1;
            uOut(k,1) = fun(vOut{end}(idx:idx+N-1,:,:), dom, 'trig');
        end
    end
else % e.g., tspan = [0, t1, t2], output the solutions at t = 0, t1 and t2
    N = length(vOut{1})/nVars;
    uOut = chebmatrix(fun(vOut{1}(1:N,:,:), dom, 'trig'));
    for k = 2:nVars
        idx = (k-1)*N + 1;
        uOut(k,1) = fun(vOut{1}(idx:idx+N-1,:,:), dom, 'trig');
    end
    for l = 2:size(vOut, 2)
        N = length(vOut{l})/nVars;
        for k = 1:nVars
            idx = (k-1)*N + 1;
            uOut(k,l) = fun(vOut{l}(idx:idx+N-1,:,:), dom, 'trig');
        end
    end
end

% Output TOUT:
if ( nargout > 2 )
    if ( length(tspan) == 2 )
        tOut = tOut(2);
    end
end

end