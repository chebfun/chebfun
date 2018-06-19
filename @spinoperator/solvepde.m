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

% SOLVEPDE has been called by SPIN/SPIN2/SPIN3/SPINSPHERE. The inputs have been 
% parsed in those files and are expeceted to be:
%
% OPTION 1.     SOLVEPDE(S, N, DT), S is a SPINOPERATOR object, N is the number
%               of grid points and DT is the time-step.
%
% OPTION 2.     SOLVEPDE(S, N, DT, PREF), PREF is a SPINPREFERENCE.

% Get the inputs:
if ( nargin == 3 ) % OPTION 1
    S = varargin{1};
    N = varargin{2};
    dt = varargin{3};
    pref = getPreference(S); % get the default preferences
elseif ( nargin == 4 ) % OPTION 2
    S = varargin{1};
    N = varargin{2};
    dt = varargin{3};
    pref = varargin{4};
end

% Dimension:
dim = getDimension(S);

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

% For PDEs on the sphere, if the constant in front of the Laplacian is real,
% use IMEX-BDF4, otherwise use LIRK4 (unless a specific scheme has been given
% by the user):
if ( isDiag(S) == 0 ) % Nondiagonal operators = operators on the sphere.
    if ( isreal(L) == 1 && isempty(pref.scheme) == 1 )
        pref.scheme = 'imexbdf4';
    elseif ( isreal(L) == 0 && isempty(pref.scheme) == 1 )
        pref.scheme = 'lirk4';
    end
end

% Create a time-stepping scheme:
schemeName = pref.scheme;
K = [];
try K = expinteg(schemeName);
catch
end
try K = imex(schemeName);
catch
end
if ( isempty(K) == 1 )
    error(['Unrecognized time-stepping scheme. See HELP/EXPINTEG and ', ...
        'HELP/IMEX for a list of available schemes.'])
end
q = K.steps;  % Number of steps of the scheme (q>1 for multistep schemes)

% Diagonal SPINOPERATOR objects (1D/2D/3D) use exponential integtrators while
% nondiagonal SPINOPERATOR objects (sphere) use IMEX schemes. Check we're using
% the right scheme:
if ( isDiag(S) == 1 ) % 1D/2D/3D
    if ( isa(K, 'expinteg') ~= 1 )
        error(['Use exponential integrators with SPIN/SPIN2/SPIN3. ', ...
            'See HELP/EXPINTEG.'])
    end
else % sphere
    if ( isa(K, 'imex') ~= 1 )
        error('Use IMEX schemes with SPINSPHERE. See HELP/IMEX.')
    end
end

% Exponential integrators use contour integrals for computing the phi-functions:
if ( isa(K, 'expinteg') == 1 )
    M = pref.M; % Number of points for complex means
else
    M = [];
end

% Set-up spatial grid for computation:
compGrid = getGrid(S, N, dom);

% Get the values to coeffs and coeffs to values transform:
v2c = getVals2CoeffsTransform(S);
c2v = getCoeffs2ValsTransform(S);

% Get the values VINIT and Fourier coeffs CINIT of the initial condition.
% Also store the nonlinear evaluation of the initial Fourier coefficients 
% in NCINIT.
if ( isDiag(S) == 1 ) % the linear part of the operartor is diagonal (1D/2D/3D)
% Note: In this case, we store the values/coeffs as a Nx1 vector in 1D,
% NxN matrix in 2D and NxNxN tensor in 3D. When there is more than one variable
% (i.e., for systems of PDEs), we store the values/coeffs as a (NVARS*N)x1 
% vector in 1D, (NVARS*N)xN matrix in 2D and (N*NVARS)xNxN tensor in 3D; NVARS
% is the number of variables. E.g., in 2D with NVARS=2:
%
%                  --------------
%                  |            |
%                  | NxN matrix |  <-- First variable
%                  |            |
%    VINIT    =    --------------
% (2*N)xN matrix   |            |
%                  | NxN matrix |  <-- Second variable
%                  |            |
%                  --------------
%
    % Initial values VINIT:
    vInit = [];
    for k = 1:nVars
        if ( dim == 1 ) % 1D
            xx = compGrid{1};
            vInit = [vInit; feval(u0{k}, xx)]; %#ok<*AGROW>
        elseif ( dim == 2 ) % 2D 
            xx = compGrid{1};
            yy = compGrid{2};
            vInit = [vInit; feval(u0{k}, xx, yy)];
        elseif ( dim == 3 ) % 3D
            xx = compGrid{1}; 
            yy = compGrid{2};
            zz = compGrid{3};
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
% Note: In this case, we store the values as above but the coeffs as a (N*N)x1 
% vector and as a (NVARS*N*N)x1 vector for systems. E.g., with NVARS=2:
%
%                   ---------
%                   |(N*N)x1|
%                   | V     |
%                   | E     |
%                   | C     |   <-- First variable
%                   | T     |
%                   | O     |
%                   | R     |
%    CINIT    =     --------- 
% (2*N*N)x1 vector  |(N*N)x1|
%                   | V     |
%                   | E     |
%                   | C     |
%                   | T     |   <-- Second variable
%                   | O     |
%                   | R     |
%                   ---------
%
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

% Get enough initial data when using a multistep scheme:
if ( q > 1 )
    [cInit, NcInit] = startMultistep(K, dt, L, Nc, Nv, pref, S, cInit, NcInit);
end

% Compute the coefficients of the time-stepping schemes:
schemeCoeffs = computeCoeffs(K, dt, L, M, S);

% Indexes for dealiasing:
ind = getDealiasingIndexes(S, N, nVars);

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
    % NPLOT is the number of grid points for the finer grid used for plotting
    % data. If NPLOT > N, create a finer grid:
    if ( pref.Nplot > N )
        plotGrid = getGrid(S, pref.Nplot, dom);
    % Otherwise, use the same grid:
    else
        plotGrid = compGrid;
    end
    % Add the (periodic) endpoints to the grid:
    compGrid = reshapeGrid(S, compGrid); 
    plotGrid = reshapeGrid(S, plotGrid); 
    [p, opts] = initializeMovie(S, dt, pref, vInit, compGrid, plotGrid);
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
    [cNew, NcNew] = oneStep(K, dt, schemeCoeffs, Nc, Nv, nVars, S, cOld, NcOld);
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
        
        % Remove small imaginary parts:
        v = [];
        for k = 1:nVars
            if ( isDiag(S) == 1 ) % 1D/2D/3D
                idx = (k-1)*N + 1;
                temp = c2v(cNew{1}(idx:idx+N-1,:,:));
            else % sphere
                idx = (k-1)*N^2 + 1;
                temp = c2v(cNew{1}(idx:idx+N^2-1));
            end
            if ( max(abs(imag(temp(:)))) < max(abs(temp(:)))*1e-10 )
                temp = real(temp);
            end
            v = [v; temp];
            valuesUpdated = 1;
        end
        
        % Check if the user gave limits for the y-axis (in 1D; PREF.YLIM)
        % or for the colorbar (in 2D/3D and on the sphere; PREF.CLIM):
        if ( isa(S, 'spinop') == 1 ) % 1D
            isLimGiven = ~isempty(pref.Ylim);
        else % 2D/3D and sphere
            isLimGiven = ~isempty(pref.Clim);
            
        end
        % If that's the case, then we're going to use these limits:
        % (See individual INITIALIZEMOVIE codes for details.)
        if ( isLimGiven == 1 )
            updateMovie(S, dt, p, opts, t, v, compGrid, plotGrid);
            
        % Otherwise, the code will automatically chose and update these limits
        % and store them in the first entry of the CELL-ARRAY OPTS:
        % (OPTS also stores other informations relative to graphics; see 
        % individual INITIALIZEMOVIE codes for details.)
        else
            opts = updateMovie(S, dt, p, opts, t, v, compGrid, plotGrid);
        end

    % Store the values every ITERPLOT iterations if using WATERFALL:
    % (Remark: Only in 1D.)
    elseif ( strcmpi(plotStyle, 'waterfall') == 1 && mod(iter, iterplot) == 0 )
        v = [];
        for k = 1:nVars
            idx = (k-1)*N + 1;
            temp = c2v(cNew{1}(idx:idx+N-1));
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
                if ( isDiag(S) == 1 ) % 1D/2D/3D
                    idx = (k-1)*N + 1;
                    temp = c2v(cNew{1}(idx:idx+N-1,:,:));
                else % sphere
                    idx = (k-1)*N^2 + 1;
                    temp = c2v(cNew{1}(idx:idx+N^2-1));
                end
                if ( max(abs(imag(temp(:)))) < max(abs(temp(:)))*1e-10 )
                    temp = real(temp);
                end
                v = [v; temp];
            end
        end
        
        % Outpute values and times:
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
    updateMovie(S, dt, p, opts, t, v, compGrid, plotGrid);
end

% Use WATERFALL if using WATERFALL:
if ( strcmpi(plotStyle, 'waterfall') == 1 ) % 1D only
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
        set(gca, 'FontSize', 12), box on
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
fun = getChebfunType(S);

% In 3D, the data comes from MESHGRID, need to permute it because the CHEBFUN3
% constructor assumes that the data comes from NDGRID:
if ( dim == 3 ) 
    for l = 1:size(vOut, 2)
        for k = 1:nVars
            idx = (k-1)*N + 1;
            vals = vOut{l}(idx:idx+N-1,:,:);
            vOut{l}(idx:idx+N-1,:,:) = permute(vals, [2 1 3]);
        end
    end
end

% Reshape the data that will be used for constructing the solution: 
% (For example, for SPINOPSPEHRE, it extracts half of the data, since the data 
% has been doubled-up with the DFS method.)
vOut = reshapeData(S, vOut, nVars);

% Output a CHEBFUN/CHEBMATRIX from values VOUT.
% Case 1; TSPAN = [0 T], output only the solution at T:
if ( length(tspan) == 2 ) 
    
    % Case 1.1; one variable:
    if ( nVars == 1 )
        uOut = fun(vOut{end}, dom, 'trig');
        
    % Case 1.2; systems, use a CHEBAMTRIX:
    else
        N = size(vOut{1}, 1)/nVars;
        uOut = chebmatrix(fun(vOut{end}(1:N,:,:), dom, 'trig'));
        for k = 2:nVars
            idx = (k-1)*N + 1;
            uOut(k,1) = fun(vOut{end}(idx:idx+N-1,:,:), dom, 'trig');
        end
    end
    
% Case 2; TSPAN = [0, t1, t2], output the solutions at t = 0, t1 and t2:
% (We store the functions in a CHEBMATRIX. The rows correspond to the different
% variables while the columns correspond to the different times.)
else
    N = size(vOut{1}, 1)/nVars;
    uOut = chebmatrix(fun(vOut{1}(1:N,:,:), dom, 'trig'));
    for k = 2:nVars
        idx = (k-1)*N + 1;
        uOut(k,1) = fun(vOut{1}(idx:idx+N-1,:,:), dom, 'trig');
    end
    for l = 2:size(vOut, 2)
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