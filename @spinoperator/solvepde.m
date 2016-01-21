function [uOut, tOut] = solvepde(varargin)
%SOLVEPDE  
%   SOLVEPDE

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%% Inputs:

% Get the inputs:
pref = [];
S = [];
tspan = [];
u0 = [];
for j = 1:nargin
    item =  varargin{j};
    if ( isa(item, 'char') == 1 )
        pdechar = item;
    elseif ( isa(item, 'double') == 1 ) 
        tspan = item;
    elseif ( isa(item, 'chebfun') == 1 )
        u0 = chebmatrix(item);
    elseif ( isa(item, 'chebfun2') == 1 )
        u0 = chebmatrix(item);
    elseif ( isa(item, 'chebfun2v') == 1 )
        u0 = chebmatrix(item(1));
        for k = 2:size(item, 1)
            u0(k,1) = item(k);
        end
    elseif ( isa(item, 'chebfun3') == 1 )
        u0 = chebmatrix(item);
    elseif ( isa(item, 'chebfun3v') == 1 )
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
end

% Parse the inputs:
if ( isempty(tspan) == 1 && isempty(u0) == 1 )
    if ( isempty(pref) == 0 )
        [tspan, u0] = parseInputs(pdechar);
    else
        [tspan, u0, pref] = parseInputs(pdechar);
    end
elseif ( isempty(tspan) == 1 )
    if ( isempty(pref) == 0 )
        tspan = parseInputs(pdechar);
    else
        [tspan, ~, pref] = parseInputs(pdechar);
    end
elseif ( isempty(u0) == 1 )
    if ( isempty(pref) == 0 )
        [~, u0] = parseInputs(pdechar);
    else
        [~, u0, pref] = parseInputs(pdechar);
    end
end

% Space interval DOM and final time TF:
dom = u0{1}.domain;
tf = tspan(end);

%% Pre-processing:

% If no SPINOPERATOR was given (i.e., DEMO mode), create one:
if ( isempty(S) == 1 ) 
    if ( isa(u0{1}, 'chebfun') == 1 )
        dim = 1;
        S = spinop(pdechar, dom);
    elseif ( isa(u0{1}, 'chebfun2') == 1 )
        dim = 2;
        S = spinop2(pdechar, dom);
    elseif ( isa(u0{1}, 'chebfun3') == 1 )
        dim = 3;
        S = spinop3(pdechar, dom);
    end
else
    dim = S.dimension;
end

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

% Dealiasing (0=NO, 1=YES):
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
    N = Nmin;
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
    dt = dtmax;
else
    % Not adpative in time, i.e., use the dt given by the user:
    dt = pref.dt;
end

% Plotting options:
plottingstyle = pref.plotting;

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
nVars = S.numVars;
xx = trigpts(N, dom(1:2));
if ( dim == 2 )
    [xx, yy] = meshgrid(xx);
elseif ( dim == 3 )
    [xx, yy, zz] = meshgrid(xx);
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
if ( strcmpi(plottingstyle, 'waterfall') == 1 )
    vWater{1} = vInit;
    twater = 0;
end

% Plot initial condition if using MOVIE:
if ( strcmpi(plottingstyle, 'movie') == 1 )
    if ( dim == 1 )
        gridpts = xx;
    elseif ( dim == 2 )
        gridpts = {xx; yy};
    elseif ( dim == 3 );
        gridpts = {xx; yy; zz};
    end
    [p, plotOptions] = initializeMovie(S, dt, pref, vInit, gridpts);
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
    if ( dealias == 1 )
        cNew{1}(ind) = 0;
    end
    
    % Check if N is large enough (i.e., check resolution in space):
    ishappy = checkHappiness(S, cNew, pref);
    
    % If resolved in space, or N>=Nmax, or not adpative in space, check if
    % resolved in time:
    if ( ishappy == 1 || N >= Nmax || adaptiveSpace == 0 )
        
        % Two steps in time with DT/2 and N points (if adpative in time):
        if ( adaptiveTime == 1 )
            [cNew2, NcNew2] = oneStep(K, schemeCoeffs2, Nc, Nv, nVars, ...
                cOld, NcOld);
            if ( dealias == 1 )
                cNew2{1}(ind) = 0;
            end
            [cNew2, NcNew2] = oneStep(K, schemeCoeffs2, Nc, Nv, nVars, ...
                cNew2, NcNew2);
            if ( dealias == 1 )
                cNew2{1}(ind) = 0;
            end
            err = max(max(max(abs(cNew{1} - cNew2{1}))));
            err = err/max(max(max(abs(cNew2{1}))));
            
        % If not adaptive in time, set CNEW2=CNEW and ERR=1:
        else
            err = 1;
            cNew2 = cNew;
            NcNew2 = NcNew;
        end
        
        % If |cnew - cnew2| is small enough, or DT<=DTMIN, or not adadptive in
        % time, keep the step:
        if (  err <= errTol || dt <= dtmin || adaptiveTime == 0 )
            
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
            if ( strcmpi(plottingstyle, 'movie') == 1 && ...
                    mod(iter,iterplot) == 0 )
                v = [];
                for k = 1:nVars
                    idx = (k-1)*N + 1;
                    temp = ifftn(cOld{1}(idx:idx+N-1,:,:));
                    if ( max(max(max(abs(imag(temp))))) < errTol )
                        temp = real(temp);
                    end
                    v = [v; temp];
                end
                valuesUpdated = 1;
                plotOptions = plotMovie(S, dt, p, plotOptions, t, v, gridpts);
                
            % Store the values every ITERPLOT iterations if using WATERFALL:
            % (Remark: Only in dimension 1.)
            elseif ( strcmpi(plottingstyle, 'waterfall') == 1 && ...
                    mod(iter, iterplot) == 0 )
                v = [];
                for k = 1:nVars
                    idx = (k-1)*N + 1;
                    temp = ifft(cOld{1}(idx:idx+N-1));
                    if ( max(max(max(abs(imag(temp))))) < errTol )
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
                        if ( max(max(max(abs(imag(temp))))) < errTol )
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
                    [xx, yy] = meshgrid(xx);
                    u = chebfun2(valsOld, dom, 'trig');   
                    temp = feval(u, xx, yy);
                elseif ( dim == 3 )
                    xx = trigpts(NN, dom(1:2));
                    [xx, yy, zz] = meshgrid(xx);
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
        
        % Update the indexes for dealiasing:
        toOne = floor(N/2) + 1 - ceil(N/6):floor(N/2) + ceil(N/6);
        if ( dim == 1 )
            ind = false(N, 1);
            ind(toOne) = 1;
            gridpts = trigpts(N, dom(1:2));
        elseif ( dim == 2 )
            ind = false(N, N);
            ind(toOne, toOne) = 1;  
            gridpts = {xx; yy};
        elseif ( dim == 3 );
            ind = false(N, N, N);
            ind(toOne, toOne, toOne) = 1;              
            gridpts = {xx; yy; zz};
        end
        ind = repmat(ind, nVars, 1);

        success = 0;
        
    end
    
end

%% Post-processing:

% Make sure that the solution at TF has been plotted if using MOVIE:
if ( strcmpi(plottingstyle, 'movie') == 1 )
    plotMovie(S, dt, p, plotOptions, t, v, gridpts);
end

% Use WATERFALL if using WATERFALL:
if ( strcmpi(plottingstyle, 'waterfall') == 1 )
    clf
    for k = 1:nVars
        uwater = [];
        for l = 1:size(vWater, 2)
            N = length(vWater{l})/nVars;
            idx = (k-1)*N + 1;
            uwater = [ uwater, chebfun(real(vWater{l}(idx:idx+N-1)), ...
                dom, 'trig') ]; %#ok<*AGROW>
        end
        subplot(nVars, 1, k)
        waterfall(uwater, twater), axis([xx(1), xx(end), 0, tf])
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

% Simpliyfy: (Remeark: CHEBFUN3/SIMPLIFY not implemented yet.)
if ( nVars == 1 )
    if ( dim == 1 || dim == 2 )
        uOut = simplify(uOut);
    end
else
    if ( dim == 1 || dim == 2 )
        doSimplify = @(f) simplify(f);
        uOut.blocks = cellfun(doSimplify, uOut.blocks, 'UniformOutput', 0);
    end
end

% Output TOUT:
if ( nargout == 2 )
    if ( length(tspan) == 2 )
        tOut = tOut(2);
    end
end

end

%% Function to parse the inputs:

function [tspan, u0, pref] = parseInputs(pdechar)
%PARSEINPUTS   Parse inputs.

if ( strcmpi(pdechar, 'AC') == 1 )
    tspan = [0 300];
    dom = [0 2*pi];
    u0 = chebmatrix(chebfun(@(x) tanh(2*sin(x)) + 3*exp(-27.*(x-4.2).^2) ...
        - 3*exp(-23.5.*(x-pi/2).^2) + 3*exp(-38.*(x-5.4).^2), dom, 'trig'));
    pref = [];
    
elseif ( strcmpi(pdechar, 'Burg') == 1 )
    tspan = [0 30];
    dom = [-1 1];
    u0 = chebmatrix(chebfun('(1-x.^2).*exp(-30.*(x+1/2).^2)', dom, 'trig'));
    pref = [];
    
elseif ( strcmpi(pdechar, 'CH') == 1 )
    tspan = [0 70];
    dom = [-1 1];
    u0 = chebmatrix(chebfun('(sin(4*pi*x)).^5 - sin(pi*x)', dom, 'trig'));
    pref = [];
    
elseif ( strcmpi(pdechar, 'GL2') == 1 )
    tspan = [0 150];
    vals = .1*randn(128, 128);
    dom = [0 200 0 200];
    u0 = chebmatrix(chebfun2(vals, dom, 'trig'));
    pref = [];
    
elseif ( strcmpi(pdechar, 'GL3') == 1 )
    tspan = [0 300];
    vals = .1*randn(32, 32, 32);
    dom = [0 100 0 100 0 100];
    u0 = chebmatrix(chebfun3(vals, dom, 'trig'));
    pref = [];
    
elseif ( strcmpi(pdechar, 'GS') == 1 )
    tspan = [0 10000];
    L = 50;
    dom = L*[-1 1];
    u01 = chebfun(@(x) 1 - 1/2*sin(pi*(x-L)/(2*L)).^100, dom, 'trig');
    u02 = chebfun(@(x) 1/4*sin(pi*(x-L)/(2*L)).^100, dom, 'trig');
    u0 = chebmatrix(u01);
    u0(2,1) = u02;
    pref = spinpref('N', 256, 'dt', 5, 'iterPlot', 10);
    
elseif ( strcmpi(pdechar, 'GS2') == 1 )
    tspan = [0 3200];
    G = 1.25;
    dom = G*[0 1 0 1];
    u01 = chebfun2(@(x,y) 1 - exp(-150*((x-G/2).^2 + (y-G/2).^2)), dom, 'trig');
    u02 = chebfun2(@(x,y) exp(-150*((x-G/2).^2 + 2*(y-G/2).^2)), dom, 'trig');
    u0 = chebmatrix(u01);
    u0(2,1) = u02;
    pref = spinpref2('dt', 8);
    
elseif ( strcmpi(pdechar, 'GS3') == 1 )
    tspan = [0 1000];
    G = 0.75;
    dom = G*[0 1 0 1 0 1];
    u01 = chebfun3(@(x,y,z) 1 - exp(-150*((x-G/2).^2 + (y-G/2).^2 + ...
        (z-G/2).^2)), dom, 'trig');
    u02 = chebfun3(@(x,y,z) exp(-150*((x-G/2).^2 + 2*(y-G/2).^2 + ...
        (z-G/2).^2)), dom, 'trig');
    u0 = chebmatrix(u01);
    u0(2,1) = u02;
    pref = spinpref3('dt', 8);
    
elseif ( strcmpi(pdechar, 'KdV') == 1 )
    A = 25^2; B = 16^2;
    tspan = [0 2*pi*3/A];
    dom = [-pi pi];
    u0 = @(x) 3*A*sech(.5*sqrt(A)*x).^2 + 3*B*sech(.5*sqrt(B)*(x-1)).^2;
    u0 = chebmatrix(chebfun(u0, dom, 'trig'));
    pref = [];
    
elseif ( strcmpi(pdechar, 'KS') == 1 )
    tspan = [0 300];
    dom = [0 32*pi];
    u0 = chebmatrix(chebfun('cos(x/16).*(1 + sin((x-1)/16))', dom, 'trig'));
    pref = [];
    
elseif ( strcmpi(pdechar, 'NLS') == 1 )
    tspan = [0 20];
    A = 2; B = 1;
    dom = [-pi pi];
    u0 = @(x) (2*B^2./(2 - sqrt(2)*sqrt(2-B^2)*cos(A*B*x)) - 1)*A;
    u0 = chebmatrix(chebfun(u0, dom, 'trig'));
    pref = [];
    
elseif ( strcmpi(pdechar, 'Schnak2') == 1 ) 
    tspan = [0 300];
    G = 50;
    dom = G*[0 1 0 1];
    u01 = chebfun2(@(x,y) 1 - exp(-10*((x-G/2).^2 + (y-G/2).^2)), dom, 'trig');
    u02 = chebfun2(@(x,y) exp(-10*((x-G/2).^2 + 2*(y-G/2).^2)), dom, 'trig');
    u0 = chebmatrix(u01);
    u0(2,1) = u02;
    pref = [];
    
elseif ( strcmpi(pdechar, 'Schnak3') == 1 )
    tspan = [0 400];
    G = 20;
    dom = G*[0 1 0 1 0 1];
    u01 = chebfun3(@(x,y,z) 1 - exp(-10*((x-G/2).^2 + (y-G/2).^2 + ...
        (z-G/2).^2)), dom, 'trig');
    u02 = chebfun3(@(x,y,z) exp(-10*((x-G/2).^2 + 2*(y-G/2).^2 + ...
        (z-G/2).^2)), dom, 'trig');
    u0 = chebmatrix(u01);
    u0(2,1) = u02;
    pref = [];
    
elseif ( strcmpi(pdechar, 'SH2') == 1 ) 
    tspan = [0 200];
    dom = [0 50 0 50];
    vals = .1*randn(64, 64);
    u0 = chebmatrix(chebfun2(vals, dom, 'trig'));
    %u0 = @(x,y) 1/4*(sin(pi*x/10) + sin(pi*y/10) + sin(pi*x/2).*sin(pi*y/2));
    %u0 = chebmatrix(chebfun2(u0, dom));
    pref = [];
    
elseif ( strcmpi(pdechar, 'SH3') == 1 ) 
    tspan = [0 200];
    dom = [0 50 0 50 0 50];
    %vals = .1*randn(32, 32, 32);
    %u0 = chebmatrix(chebfun3(vals, dom, 'trig'));
    u0 = @(x,y,z) 1/4*(sin(pi*x/10) + sin(pi*y/10) + sin(pi*z/10) + ...
        sin(pi*x/2).*sin(pi*y/2) + sin(pi*x/2).*sin(pi*z/2) + ...
        sin(pi*z/2).*sin(pi*y/2));
    u0 = chebmatrix(chebfun3(u0, dom));
    pref = [];
    
else
    error('SPINOPERATOR:SOLVEPDE:parseInputs', 'Unrecognized PDE.')
end

end