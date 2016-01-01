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

%%
tout = [];
uout = [];

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