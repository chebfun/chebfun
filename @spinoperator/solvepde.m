function [uOut, tOut] = solvepde(varargin)
%SOLVEPDE   Solve a PDE defined by a SPINOP, a SPINOP2 or a SPINOP3.
%   SOLVEPDE is called by SPIN, SPIN2 and SPIN3. It is not called directly by
%   the user. Appropriate help texts can be found in SPIN, SPIN2 and SPIN3.
%
% See also SPIN, SPIN2, SPIN3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% Inputs:

% Throw an error if no input:
if ( nargin == 0 ) 
    error('SPINOPERATOR:solvepde', 'Not enough input arguments.')
end

% Throw an error if the first input is not appropriate:
if ( nargin == 1 )
   item = varargin{1};
   if ( isa(item, 'spinoperator') == 0 && isa(item, 'char') == 0 )
       error('SPINOPERATOR:solvepde', ['Firt input should be a ', ...
           'SPINOPERATOR or a STRING.'])
   end
end

% Get the inputs:
pref = [];
S = [];
tspan = [];
u0 = [];
j = 1;
while ( j <= nargin )
    item =  varargin{j};
    if ( isa(item, 'char') == 1 )
        pdechar = item;
    elseif ( isa(item, 'double') == 1 ) 
        tspan = item;
    elseif ( isa(item, 'chebfun') == 1 || isa(item, 'chebfun2') == 1 || ...
        isa(item, 'chebfun3') == 1 )
        u0 = chebmatrix(item);
    elseif ( isa(item, 'chebfun2v') == 1 || isa(item, 'chebfun3v') == 1 )
        u0 = chebmatrix(item(1));
        for k = 2:size(item, 1)
            u0(k,1) = item(k);
        end
    elseif ( isa(item, 'chebmatrix') == 1 )
        u0 = item;
    elseif ( isa(item, 'spinoperator') == 1 )
        S = item;
    elseif ( isa(item, 'spinpreference') == 1 )
        pref = item;
    else
        error('SPINOPERATOR:solvepde', 'Unrecognized input.')
    end
    j = j + 1;
end

% A SPINOPERATOR was given by the user:
if ( isempty(S) == 0 )
    dim = getDimension(S);
    
    % Create a SPINPREFERENCE object if none:
    if ( isempty(pref) == 1 )
        if ( dim == 1 )
            pref = spinpref();
        elseif ( dim == 2 )
            pref = spinpref2();
        elseif ( dim == 3 )
            pref = spinpref3();
        end
    end
    
% DEMO mode, i.e., a STRING was given by the user:    
else
    
    % Create a SPINOPERATOR:
    is2D = ~isempty(strfind(pdechar, '2'));
    is3D = ~isempty(strfind(pdechar, '3'));
    is1D = ( is2D == 0 && is3D == 0 );
    if ( is1D == 1 )
        dim = 1;
        S = spinop(pdechar);
    elseif ( is2D == 1 )
        dim = 2;
        S = spinop2(pdechar);
    elseif ( is3D == 1 )
        dim = 3;
        S = spinop3(pdechar);
    end
    
    % Create a SPINPREFERENCE object if none:
    if ( isempty(pref) == 1 )
        if ( dim == 1 )
            if ( isempty(tspan) == 1 )
                pref = spinpref(pdechar);
            else
                pref = spinpref();
            end
        elseif ( dim == 2 )
            if ( isempty(tspan) == 1 )
                pref = spinpref2(pdechar);
            else
                pref = spinpref2();
            end
        elseif ( dim == 3 )
            if ( isempty(tspan) == 1 )
                pref = spinpref3(pdechar);
            else
                pref = spinpref3();
            end
        end
    end
end
   
% Time interval TSPAN:
if ( isempty(tspan) == 1 )
    tspan = S.tspan;
end
    
% Initial condition U0:
if ( isempty(u0) == 1 )
    u0 = S.init;
    if ( isa(u0, 'chebfun') == 1 || isa(u0, 'chebfun2') == 1 || ...
        isa(u0, 'chebfun3') == 1 )
        u0 = chebmatrix(u0);
    elseif ( isa(u0, 'chebfun2v') == 1 || isa(u0, 'chebfun3v') == 1 )
        temp = chebmatrix(u0(1));
        for k = 2:size(u0, 1)
            temp(k,1) = u0(k);
        end
        u0 = temp;
    end
end

% Convert to trigfun:
nVars = S.numVars;
for k = 1:nVars
    if ( dim == 1 )
        temp = chebfun(u0{k}, 'trig');
    elseif ( dim == 2 )
        temp = chebfun2(u0{k}, 'trig');
    elseif ( dim == 3 )
        temp = chebfun3(u0{k}, 'trig');
    end
    u0(k,1) = temp;
end

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

%% Pre-processing:

% Dealiasing:
dealias = pref.dealias;

% Error tolerance:
errTol = pref.errTol;

% Points for complex means:
M = pref.M;

% Plot every ITERPLOT iterations if 'movie':
iterplot = pref.iterPlot;

% Number of grid points N:
Nmin = pref.Nmin;
Nmax = pref.Nmax;
adaptiveSpace = isempty(pref.N);
if ( adaptiveSpace == 1 )
    % Adaptive in space, start with NMIN:
    % (Unless NMIN is smaller than the length NU0 of the initial condition.)
    Nu0 = max(cellfun(@(B) length(B), u0.blocks)); % Length
    Nu0 = 2^ceil(log2(Nu0)); % Convert to a power of 2
    Nu0 = min(Nu0, Nmax); % Make sure it's not larger than Nmax
    N = max(Nmin, Nu0);
else
    % Not adpative in space, i.e., use the N given by the user:
    N = pref.N;
end

% Time-step dt:
dtmin = pref.dtmin;
dtmax = pref.dtmax;
adaptiveTime = isempty(pref.dt);
if ( adaptiveTime == 1 )
    % Adaptive in time, start with DTMAX: 
    % (Unless the interval is shorter than DTMAX.)
    dt = min(tspan(end), dtmax);
else
    % Not adpative in time, i.e., use the dt given by the user:
    dt = pref.dt;
end

% Plot options:
plotStyle = pref.plot;

% Create a time-stepping scheme:
schemeName = pref.scheme;
K = spinscheme(schemeName);

% Get the number of steps of the scheme:
q = K.steps;

% Operators (linear part L, and nonlinear parts Nc and Nv):
[L, Nc] = discretize(S, N);
Nv = S.nonlinearPartVals;

% Set-up spatial grid, and initial condition (values VINIT and Fourier coeffs 
% CINIT):
xx = trigpts(N, dom(1:2));
if ( dim == 2 )
    yy = trigpts(N, dom(3:4));
    [xx, yy] = meshgrid(xx, yy);
elseif ( dim == 3 )
    yy = trigpts(N, dom(3:4));
    zz = trigpts(N, dom(5:6));
    [xx, yy, zz] = meshgrid(xx, yy, zz);
end
vInit = [];
for k = 1:nVars
    if ( dim == 1 )
        vInit = [vInit; feval(u0{k}, xx)]; %#ok<*AGROW>
    elseif ( dim == 2 )
        vInit = [vInit; feval(u0{k}, xx, yy)];
    elseif ( dim == 3 )
        vInit = [vInit; feval(u0{k}, xx, yy, zz)];
    end
end
cInit{1} = [];
for k = 1:nVars
    idx = (k-1)*N + 1;
    cInit{1} = [cInit{1}; fftn(vInit(idx:idx+N-1,:,:))];
end

% Store the nonlinear evaluation of the initial data in NCINIT:
vals = ifftn(cInit{1}(1:N,:,:));
for k = 1:nVars-1
    idx = k*N + 1;
    vals = [vals; ifftn(cInit{1}(idx:idx+N-1,:,:))];
end
vals = Nv(vals);
coeffs = fftn(vals(1:N,:,:));
for k = 1:nVars-1
    idx = k*N + 1;
    coeffs = [coeffs; fftn(vals(idx:idx+N-1,:,:))];
end
coeffs = Nc.*coeffs;
NcInit{1} = coeffs;
    
% Get enough initial data when using a multistep scheme:
if ( q > 1 )
    [cInit, NcInit, dt] = startMultistep(K, adaptiveTime, dt, L, Nc, Nv, ...
        pref, S, cInit, NcInit);
end
vInit = ifftn(cInit{1}(1:N,:,:));
for k = 1:nVars-1
    idx = k*N + 1;
    vInit = [vInit; ifftn(cInit{1}(idx:idx+N-1,:,:))];
end

% Compute the coefficients of the scheme:
schemeCoeffs = computeCoeffs(K, dt, L, M, S);

% If adaptive in time, get the coefficients with DT/2:
if ( adaptiveTime == 1 )
    schemeCoeffs2 = computeCoeffs(K, dt/2, L, M, S);
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
        dataGrid = {xx};
        plotGrid = {trigpts(Nplot, dom)};
    elseif ( dim == 2 )
        dataGrid = {xx; yy};
        ttx = trigpts(Nplot, dom(1:2));
        tty = trigpts(Nplot, dom(3:4));
        [xxx, yyy] = meshgrid(ttx, tty);
        plotGrid = {xxx; yyy};
    elseif ( dim == 3 );
        dataGrid = {xx; yy; zz};
        ttx = trigpts(Nplot, dom(1:2));
        tty = trigpts(Nplot, dom(3:4));
        ttz = trigpts(Nplot, dom(5:6));
        [xxx, yyy, zzz] = meshgrid(ttx, tty, ttz);
        plotGrid = {xxx; yyy; zzz};
    end
    dataGrid = reshapeGrid(S, dataGrid);
    plotGrid = reshapeGrid(S, plotGrid);
    [p, options] = initializeMovie(S, dt, pref, vInit, dataGrid, plotGrid);
end

%% Time-stepping loop:

iter = 0;
t = (q-1)*dt;
success = 0;
pos = 2;
cOld = cInit;
NcOld = NcInit;
while ( t < tf )

    % One step in time with DT and N points:
    [cNew, NcNew] = oneStep(K, schemeCoeffs, Nc, Nv, nVars, cOld, NcOld);
    valuesUpdated = 0;
    
    % Dealiasing procedure:
    if ( strcmpi(dealias, 'on') == 1 )
        cNew{1}(ind) = 0;
    end
    
    % Check if N is large enough (i.e., check resolution in space):
    if ( adaptiveSpace == 1 )
        ishappySpace = checkHappiness(S, cNew, pref);
        
    % If not adaptive in space, set ISHAPPYSPACE=1:
    else
        ishappySpace = 1;
    end
    
    % If happy with the resolution in space or N>=Nmax:
    if ( ishappySpace == 1 || N >= Nmax )
        
        % Two steps in time with DT/2 and N points (if adpative in time):
        if ( adaptiveTime == 1 )
            [cNew2, NcNew2] = oneStep(K, schemeCoeffs2, Nc, Nv, nVars, ...
                cOld, NcOld);
            if ( strcmpi(dealias, 'on') == 1 )
                cNew2{1}(ind) = 0;
            end
            [cNew2, NcNew2] = oneStep(K, schemeCoeffs2, Nc, Nv, nVars, ...
                cNew2, NcNew2);
            if ( strcmpi(dealias, 'on') == 1 )
                cNew2{1}(ind) = 0;
            end
            err = max(abs(cNew{1}(:) - cNew2{1}(:)));
            err = err/max(abs(cNew2{1}(:)));
            ishappyTime = ( err <= errTol );
            
        % If not adaptive in time, set CNEW2=CNEW and ISHAPPYTIME=1:
        else
            cNew2 = cNew;
            NcNew2 = NcNew;
            ishappyTime = 1;
        end
        
        % If happy with the resolution in time or DT<=DTMIN:
        if (  ishappyTime == 1 || dt <= dtmin )
            
            % Update time T, iteration ITER and Fourier coefficients C:
            iter = iter + 1;
            if ( adaptiveTime == 1 )
                t = t + dt;
            else
                t = (iter + q - 1)*dt;
            end
            success = success + 1;
            cOld = cNew2;
            NcOld = NcNew2;
            
            % Plot every ITERPLOT iterations if using MOVIE:
            if ( strcmpi(plotStyle, 'movie') == 1 && ...
                    mod(iter,iterplot) == 0 )
                v = [];
                for k = 1:nVars
                    idx = (k-1)*N + 1;
                    temp = ifftn(cOld{1}(idx:idx+N-1,:,:));
                    if ( max(abs(imag(temp(:)))) < errTol )
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
                    plotMovie(S, dt, p, options, t, v, dataGrid, plotGrid);
                else
                    options = plotMovie(S, dt, p, options, t, v, dataGrid, ...
                        plotGrid);
                end
                
            % Store the values every ITERPLOT iterations if using WATERFALL:
            % (Remark: Only in dimension 1.)
            elseif ( strcmpi(plotStyle, 'waterfall') == 1 && ...
                    mod(iter, iterplot) == 0 )
                v = [];
                for k = 1:nVars
                    idx = (k-1)*N + 1;
                    temp = ifft(cOld{1}(idx:idx+N-1));
                    if ( max(abs(imag(temp(:)))) < errTol )
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
                        temp = ifftn(cOld{1}(idx:idx+N-1,:,:));
                        if ( max(abs(imag(temp(:)))) < errTol )
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
          
            % Make sure that the solution is computed at the entries of TSPAN:
            if ( t + 2*dt > tspan(pos) && t ~= tspan(pos) )
                if ( abs(t + dt - tspan(pos)) < 1e-10 || adaptiveTime == 0 )
                    continue
                else
                    dt = (tspan(pos) - t)/2;
                    schemeCoeffs = computeCoeffs(K, dt, L, M, S);
                    schemeCoeffs2 = computeCoeffs(K, dt/2, L, M, S);
                    success = 0;
                    continue
                end
            end

            % If 50 consecutive steps have been successful, double DT, and
            % update quantities which depend on DT:
            if ( adaptiveTime == 1 && success >= 50 && 2*dt < dtmax )
                dt = 2*dt;
                schemeCoeffs2 = schemeCoeffs;
                schemeCoeffs = computeCoeffs(K, dt, L, M, S);
                success = 0;
            end
            
        % If |cnew - cnew2| is not small enough, half DT, update quantities
        % which depend on DT, and redo the step:
        else
            
            % Half DT only if DT/2>DTMIN:
            if ( dt/2 > dtmin )
                dt = dt/2;
                
            % Othewise use DTMIN:
            else
                dt = dtmin;
                warning(['The time-step has been set to its minimum ', ...
                    'value %1.1e.\n'], dtmin)
            end
            
            % Update quantities:
            schemeCoeffs = schemeCoeffs2;
            schemeCoeffs2 = computeCoeffs(K, dt/2, L, M, S);
            success = 0;
            
        end
        
    % If not resolved in space, increase N, update quantities which depend on N,
    % and redo the step:
    else
        
        % Increase N:
        if ( dim == 1 )
            NN = 2*N;
        elseif ( dim == 2 )
            NN = 1.5*N;
        elseif ( dim == 3 )
            NN = 1.25*N;
        end
        
        % Update L and NC:
        [L, Nc] = discretize(S, NN);
        
        % Loop over the number of steps of the method:
        for i = 1:q
            coeffs = [];
            vals = [];
            for k = 1:nVars
                idx = (k-1)*N + 1;
                valsOld = ifftn(cOld{i}(idx:idx+N-1,:,:));
                if ( dim == 1 )
                    xx = trigpts(NN);
                    u = trigtech({valsOld, trigtech.vals2coeffs(valsOld)});
                    temp = feval(u, xx);
                elseif ( dim == 2 )
                    xx = trigpts(NN, dom(1:2));
                    yy = trigpts(NN, dom(3:4));
                    [xx, yy] = meshgrid(xx, yy);
                    u = chebfun2(valsOld, dom, 'trig');   
                    temp = feval(u, xx, yy);
                elseif ( dim == 3 )
                    xx = trigpts(NN, dom(1:2));
                    yy = trigpts(NN, dom(3:4));
                    zz = trigpts(NN, dom(5:6));
                    [xx, yy, zz] = meshgrid(xx, yy, zz);
                    u = chebfun3(valsOld, dom, 'trig');
                    temp = feval(u, xx, yy, zz);
                end
                vals = [vals; temp];
                coeffs = [coeffs; fftn(temp)];
            end
            
            % Update the Fourier coefficients:
            cOld{i} = coeffs;
            vals = Nv(vals);
            coeffs = fftn(vals(1:NN,:,:));
            for k = 1:nVars-1
                idx = k*NN + 1;
                coeffs = [coeffs; fftn(vals(idx:idx+NN-1,:,:))];
            end
            
            % Update the nonlinear evaluations:
            NcOld{i} = Nc.*coeffs;
            
        end
        
        % Compute the new coefficients for the scheme:
        N = NN;
        schemeCoeffs = computeCoeffs(K, dt, L, M, S);
        schemeCoeffs2 = computeCoeffs(K, dt/2, L, M, S);
        
        % Update the grid for plotting and the indexes for dealiasing:
        toOne = floor(N/2) + 1 - ceil(N/6):floor(N/2) + ceil(N/6);
        if ( dim == 1 )
            ind = false(N, 1);
            ind(toOne) = 1;
            dataGrid = {trigpts(N, dom(1:2))};
        elseif ( dim == 2 )
            ind = false(N, N);
            ind(toOne, toOne) = 1;  
            dataGrid = {xx; yy};
        elseif ( dim == 3 );
            ind = false(N, N, N);
            ind(toOne, toOne, toOne) = 1;              
            dataGrid = {xx; yy; zz};
        end
        dataGrid = reshapeGrid(S, dataGrid);
        ind = repmat(ind, nVars, 1);
        success = 0;
        
    end
    
end

%% Post-processing:

% Make sure that the solution at TF has been plotted if using MOVIE:
if ( strcmpi(plotStyle, 'movie') == 1 )
    plotMovie(S, dt, p, options, t, v, dataGrid, plotGrid);
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
if ( length(tspan) == 2 )
    N = length(vOut{end})/nVars;
    if ( nVars == 1 )
        uOut = fun(vOut{end}(1:N,:,:), dom, 'trig');
    else
        uOut = chebmatrix(fun(vOut{end}(1:N,:,:), dom, 'trig'));
        for k = 2:nVars
            idx = (k-1)*N + 1;
            uOut(k,1) = fun(vOut{end}(idx:idx+N-1,:,:), dom, 'trig');
        end
    end
else
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

% Simpliyfy:
if ( nVars == 1 )
    uOut = simplify(uOut);
else
    doSimplify = @(f) simplify(f, errTol);
    uOut.blocks = cellfun(doSimplify, uOut.blocks, 'UniformOutput', 0);
end

% Output TOUT:
if ( nargout > 2 )
    if ( length(tspan) == 2 )
        tOut = tOut(2);
    end
end

end