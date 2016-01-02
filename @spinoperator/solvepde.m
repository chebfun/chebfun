function [uout, tout] = solvepde(varargin)
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
varargin = varargin{1};
nargin = length(varargin);
for j = 1:nargin
    item =  varargin{j};
    if ( isa(item, 'char') == 1 )
        pdechar = item;
    elseif ( isa(item, 'double') ) 
        tspan = item;
    elseif ( isa(item, 'chebfun') )
        u0 = chebmatrix(item);
    elseif ( isa(item, 'chebfun2') )
        u0 = chebmatrix(item);
    elseif ( isa(item, 'chebfun2v') )
        u0 = chebmatrix(item(1));
        for k = 2:size(item, 1)
            u0(k,1) = item(k);
        end
    elseif ( isa(item, 'chebfun3') )
        u0 = chebmatrix(item);
    elseif ( isa(item, 'chebmatrix') )
        u0 = item;
    elseif ( isa(item, 'spinoperator') )
        S = item;
    elseif ( isa(item, 'spinpreference') )
        pref = item;
    else
        error('SPINOPERATOR:solvepde', 'Unrecognized input.')
    end
end

% Parse the inputs:
if ( isempty(tspan) == 1 && isempty(u0) == 1 )
    [tspan, u0] = parseInputs(pdechar);
elseif ( isempty(tspan) == 1 ) 
    tspan = parseInputs(pdechar);
elseif ( isempty(u0) == 1 )
    [~, u0] = parseInputs(pdechar);
end

%% Pre-processing:

% Space interval DOM and final time TF:
dom = u0{1}.domain;
tf = tspan(end);

% Create a SPINOPERATOR if none, and get the dimension DIM:
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

% Operators (linear part L, differentiation term NC of the nonliner part):
[L, Nc] = discretize(S, N);

% Compute coefficients for the time-stepping scheme:
LR = computeLR(S, dt, L, M, N);
schemeCoeffs = computeCoeffs(K, dt, L, LR, S);

% If adaptive in time, get the coefficients with DT/2:
if ( adaptiveTime == 1 )
    LR2 = computeLR(S, dt/2, L, M, N);
    schemeCoeffs2 = computeCoeffs(K, dt/2, L, LR2, S);
end

% Set-up spatial grid, and initial condition:
nVars = S.numVars;
xx = trigpts(N, dom(1:2));
if ( dim == 2 )
    [xx, yy] = meshgrid(xx, xx);
elseif ( dim == 3 )
    [xx, yy, zz] = meshgrid(xx, xx, xx);
end
v = [];
for k = 1:nVars
    if ( dim == 1 )
        v = [v; feval(u0{k}, xx)]; %#ok<*AGROW>
    elseif ( dim == 2 )
        v = [v; feval(u0{k}, xx, yy)];
    elseif ( dim == 3 )
        v = [v; feval(u0{k}, xx, yy, zz)];
    end
end
c{1} = [];
for k = 1:nVars
    idx = (k-1)*N + 1;
    c{1} = [c{1}; fftn(v(idx:idx+N-1,:,:))];
end

% Get enough initial data when using a multistep scheme:
if ( q > 1 )
    [c, dt] = startMultistep(K, adaptiveTime, dt, L, LR, Nc, pref, S, c);
end

% Plot initial condition if using MOVIE:
if ( strcmpi(plottingstyle, 'movie') == 1 )
    if ( dim == 1 )
        gridPoints = xx;
    elseif ( dim == 2 )
        gridPoints = {xx; yy};
    elseif ( dim == 3 );
        gridPoints = {xx; yy; zz};
    end
    [p, plotOption] = initializeMovie(S, dt, pref, v, gridPoints);
end


% Indexes for dealiasing:
toOne = floor(N/2)+1-ceil(N/6):floor(N/2)+ceil(N/6);
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

% Values to output:
vout{1} = v;

% Values to plot if using WATERFALL:
if ( strcmpi(plottingstyle, 'waterfall') == 1 )
    vwater{1} = v;
    twater = 0;
end

%% Time-stepping loop:

iter = 0;
t = (q-1)*dt;
success = 0;
pos = 2;
while ( t < tf )

    % One step in time with DT and N points:
    cnew = oneStep(K, schemeCoeffs, Nc, S, c);
    valuesUpdated = 0;
    
    % Dealiasing procedure:
    if ( dealias == 1 )
        cnew{1}(ind) = 0;
    end
    
    % Check if N is large enough (i.e., check resolution in space):
%     for k = 1:nVars
%         idx = (k-1)*N + 1;
%         cnewShift = ifftshift(cnew{1}(idx:idx+N-1));
%         if ( mod(N,2) == 0 )
%             cnewShift = [cnewShift(N); cnewShift(N-1:-1:N/2+1) + ...
%                 cnewShift(1:N/2-1); cnewShift(N/2)];
%         else
%             cnewShift = [cnewShift(N:-1:(N+1)/2+1,:) + cnewShift(1:(N+1)/2-1,:); ...
%                 cnewShift((N+1)/2,:)];
%         end
%         cnewShift = flipud(cnewShift);
%         cnewShift = [cnewShift(1,1); kron(cnewShift(2:end,1), [1; 1])];
%         cutoff = standardChop(cnewShift, errTol);
%         ishappy(k) = ( cutoff < N );
%     end
%     ishappy = all(ishappy == 1);
    ishappy = 1;
    
    % If resolved in space, or N>=Nmax, or not adpative in space, check if
    % resolved in time:
    if ( ishappy == 1 || N >= Nmax || adaptiveSpace == 0 )
        
        % Two steps in time with DT/2 and N points (if adpative in time):
        if ( adaptiveTime == 1 )
            cnew2 = oneStep(K, schemeCoeffs2, Nc, S, c);
            if ( dealias == 1 )
                cnew2{1}(ind) = 0;
            end
            cnew2 = oneStep(K, schemeCoeffs2, Nc, S, cnew2);
            if ( dealias == 1 )
                cnew2{1}(ind) = 0;
            end
            err = max(abs(cnew{1} - cnew2{1}))/max(abs(cnew2{1}));
            
        % If not adaptive in time, set CNEW2=CNEW and ERR=1:
        else
            err = 1;
            cnew2 = cnew;
        end
        
        % If |cnew - cnew2| is small enough, or DT<=DTMIN, or not adadptive in
        % time, keep the step:
        if (  err <= errTol || dt <= dtmin || adaptiveTime == 0 )
            
            % Update time T, iteration ITER and Fourier coefficients C:
            t = t + dt;
            iter = iter + 1;
            success = success + 1;
            c = cnew2;
            
            % Plot every ITERPLOT iterations if using MOVIE:
            if ( strcmpi(plottingstyle, 'movie') == 1 && ...
                    mod(iter,iterplot) == 0 )
                v = [];
                for k = 1:nVars
                    idx = (k-1)*N + 1;
                    temp = ifftn(c{1}(idx:idx+N-1,:,:));
                    if ( norm(imag(temp), inf) < errTol )
                        temp = real(temp);
                    end
                    v = [v; temp];
                end
                valuesUpdated = 1;
                plotOption = plotMovie(S, dt, p, plotOption, t, v, gridPoints);
                
            % Store the values every ITERPLOT iterations if using WATERFALL:
            % Note: Only in dimension 1.
            elseif ( strcmpi(plottingstyle, 'waterfall') == 1 && ...
                    mod(iter,iterplot) == 0 )
                v = [];
                for k = 1:nVars
                    idx = (k-1)*N + 1;
                    temp = ifft(c{1}(idx:idx+N-1));
                    if ( norm(imag(temp), inf) < errTol )
                        temp = real(temp);
                    end
                    v = [v; temp];
                end
                valuesUpdated = 1;
                vwater{iter/iterplot + 1} = v;
                twater = [twater, t];
            end
            
            % Make sure that the solution is computed at the entries of TSPAN:
            if ( t + 2*dt > tspan(pos) && t ~= tspan(pos) )
                dt = tspan(pos) - t;
                LR = computeLR(S, dt, L, M, N);
                schemeCoeffs = computeCoeffs(K, dt, L, LR, S);
                LR2 = computeLR(S, dt/2, L, M, N);
                schemeCoeffs2 = computeCoeffs(K, dt/2, L, LR2, S);
                success = 0;
                continue
            end
            
            % Output the solution if T correponds to an entry of TSPAN:
            if ( t == tspan(pos) )
                if ( valuesUpdated == 0 )
                    v = [];
                    for k = 1:nVars
                        idx = (k-1)*N + 1;
                        temp = ifftn(c{1}(idx:idx+N-1,:,:));
                        if ( norm(imag(temp), inf) < errTol )
                            temp = real(temp);
                        end
                        v = [v; temp];
                    end
                end
                vout{pos} = v;
                pos = pos + 1;
            end
            
            % If 50 consecutive steps have been successful, double DT, and
            % update quantities which depend on DT:
            if ( adaptiveTime == 1 && success == 50 && 2*dt < dtmax )
                dt = 2*dt;
                schemeCoeffs2 = schemeCoeffs;
                LR = computeLR(S, dt, L, M, N);
                schemeCoeffs = computeCoeffs(K, dt, L, LR, S);
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
            LR2 = computeLR(S, dt/2, L, M, N);
            schemeCoeffs2 = computeCoeffs(K, dt/2, L, LR2, S);
            success = 0;
        end
        
    % If not resolved in space, double N, update quantities which depend on N,
    % and redo the step:
    else
        for i = 1:q
            coeffs = [];
            for k = 1:nVars
                idx = (k-1)*N + 1;
                vals = ifft(c{i}(idx:idx+N-1));
                u = trigtech({vals, trigtech.vals2coeffs(vals)});
                coeffs = [coeffs; fft(feval(u, trigpts(2*N)))];
            end
            c{i} = coeffs;
        end
        N = 2*N;
        [L, Nc] = discretize(S, N);
        LR = computeLR(S, dt, L, M, N);
        schemeCoeffs = computeCoeffs(K, dt, L, LR, S);
        LR2 = computeLR(S, dt/2, L, M, N);
        schemeCoeffs2 = computeCoeffs(K, dt/2, L, LR2, S);
        ind = false(N, 1);
        ind(floor(N/2)+1-ceil(N/6):floor(N/2)+ceil(N/6)) = 1;
        xx = trigpts(N, dom);
        if ( dim == 1 )
            gridPoints = xx;
        elseif ( dim == 2 )
            gridPoints = {xx; yy};
        elseif ( dim == 3 );
            gridPoints = {xx; yy; zz};
        end
        success = 0;
    end
    
end

% Make sure that the solution at TF has been plotted if using MOVIE:
if ( strcmpi(plottingstyle, 'movie') == 1 )
    plotMovie(S, dt, p, plotOption, t, v, gridPoints);
end

% Use WATERFALL if using WATERFALL:
if ( strcmpi(plottingstyle, 'waterfall') == 1 )
    clf
    for k = 1:nVars
        uwater = [];
        for l = 1:size(vwater, 2)
            N = length(vwater{l})/nVars;
            idx = (k-1)*N + 1;
            uwater = [ uwater, chebfun(real(vwater{l}(idx:idx+N-1)), ...
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

%% Post-processing:

% Get the right type of CHEBFUN:
if ( dim == 1 )
    fun = @chebfun;
elseif ( dim == 2 )
    fun = @chebfun2;
elseif ( dim == 3 )
    fun = @chebfun3;
end

% Output a CHEBFUN/CHEBMATRIX from values VOUT:
if ( length(tspan) == 2 )
    N = length(vout{end})/nVars;
    if ( nVars == 1 )
        uout = fun(vout{end}(1:N,:,:), dom, 'trig');
    else
        uout = chebmatrix(fun(vout{end}(1:N,:,:), dom, 'trig'));
        for k = 2:nVars
            idx = (k-1)*N + 1;
            uout(k,1) = fun(vout{end}(idx:idx+N-1,:,:), dom, 'trig');
        end
    end
else
    N = length(vout{1})/nVars;
    uout = chebmatrix(fun(vout{1}(1:N,:,:), dom, 'trig')); 
    for k = 2:nVars
        idx = (k-1)*N + 1;
        uout(k,1) = fun(vout{1}(idx:idx+N-1,:,:), dom, 'trig');
    end
    for l = 2:size(vout, 2)
        N = length(vout{l})/nVars;
        for k = 1:nVars
            idx = (k-1)*N + 1;
            uout(k,l) = fun(vout{l}(idx:idx+N-1,:,:), dom, 'trig'); 
        end
    end
end

% Output TOUT:
if ( nargout == 2 )
    if ( length(tspan) == 2 )
        tout = tspan(2);
    else
        tout = tspan;
    end
end

end

%% Function to parse the inputs:

function [tspan, u0] = parseInputs(pdechar)
%PARSEINPUTS   Parse inputs.

if ( strcmpi(pdechar, 'AC') == 1 )
    tspan = [0 300];
    dom = [0 2*pi];
    u0 = chebmatrix(chebfun(@(x) tanh(2*sin(x)) + 3*exp(-27.*(x-4.2).^2) ...
        - 3*exp(-23.5.*(x-pi/2).^2) + 3*exp(-38.*(x-5.4).^2), dom, 'trig'));
    
elseif ( strcmpi(pdechar, 'Burg') == 1 )
    tspan = [0 30];
    dom = [-1 1];
    u0 = chebmatrix(chebfun('(1-x.^2).*exp(-30.*(x+1/2).^2)', dom, 'trig'));

elseif ( strcmpi(pdechar, 'CH') == 1 )
    tspan = [0 70];
    dom = [-1 1];
    u0 = chebmatrix(chebfun('(sin(4*pi*x)).^5 - sin(pi*x)', dom, 'trig'));

elseif ( strcmpi(pdechar, 'GL2') == 1 )
    tspan = [0 150];
    vals = .1*randn(128, 128);
    dom = [0 200 0 200];
    u0 = chebmatrix(chebfun2(vals, dom, 'trig'));
    
elseif ( strcmpi(pdechar, 'GL3') == 1 )
    tspan = [0 150];
    vals = .1*randn(32, 32, 32);
    dom = [0 100 0 100 0 100];
    u0 = chebmatrix(chebfun3(vals, dom, 'trig'));

elseif ( strcmpi(pdechar, 'GS') == 1 )
    tspan = [0 2000];
    L = 50;
    dom = L*[-1 1];
    u01 = chebfun(@(x) 1 - 1/2*sin(pi*(x-L)/(2*L)).^100, dom, 'trig');
    u02 = chebfun(@(x) 1/4*sin(pi*(x-L)/(2*L)).^100, dom, 'trig');
    u0 = chebmatrix(u01);
    u0(2,1) = u02;

elseif ( strcmpi(pdechar, 'GS2') == 1 )
    tspan = [0 20000];
    G = 1.25;
    dom = G*[0 1 0 1];
    u01 = chebfun2(@(x,y) 1 - exp(-150*((x-G/2).^2 + (y-G/2).^2)), dom, 'trig');
    u02 = chebfun2(@(x,y) exp(-150*((x-G/2).^2 + 2*(y-G/2).^2)), dom, 'trig');
    u0 = chebmatrix(u01);
    u0(2,1) = u02;
    
elseif ( strcmpi(pdechar, 'GS3') == 1 )
    tspan = [0 150];
    vals = .1*randn(32, 32, 32);
    dom = [0 100 0 100 0 100];
    u01 = chebmatrix(chebfun3(vals, dom, 'trig'));
    u02 = u01;
    u0 = chebmatrix(u01);
    u0(2,1) = u02;

elseif ( strcmpi(pdechar, 'KdV') == 1 )
    A = 25^2; B = 16^2;
    tspan = [0 2*pi*3/A];
    dom = [-pi pi];
    u0 = @(x) 3*A*sech(.5*sqrt(A)*x).^2 + 3*B*sech(.5*sqrt(B)*(x-1)).^2;
    u0 = chebmatrix(chebfun(u0, dom, 'trig'));
 
elseif ( strcmpi(pdechar, 'KS') == 1 )
    tspan = [0 300];
    dom = [0 32*pi];
    u0 = chebmatrix(chebfun('cos(x/16).*(1 + sin((x-1)/16))', dom, 'trig'));
    
elseif ( strcmpi(pdechar, 'NLS') == 1 )
    tspan = [0 20];
    A = 2; B = 1;
    dom = [-pi pi];
    u0 = @(x) (2*B^2./(2 - sqrt(2)*sqrt(2-B^2)*cos(A*B*x)) - 1)*A;
    u0 = chebmatrix(chebfun(u0, dom, 'trig'));
    
elseif ( strcmpi(pdechar, 'Schnak2') == 1 ) 
    tspan = [0 1000];
    G = 50;
    dom = G*[0 1 0 1];
    u01 = chebfun2(@(x,y) 1 - exp(-10*((x-G/2).^2 + (y-G/2).^2)), dom, 'trig');
    u02 = chebfun2(@(x,y) exp(-10*((x-G/2).^2 + 2*(y-G/2).^2)), dom, 'trig');
    u0 = chebmatrix(u01);
    u0(2,1) = u02;
    
else
    error('SPINOPERATOR:SOLVEPDE:parseInputs', 'Unrecognized PDE.')
end

end