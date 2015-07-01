function pass = test_determineDiscretization(pref)
% Test DETERMINEDISCRETIZATION method:

if ( nargin == 0 )
    pref = cheboppref();
end

%% Test the VALUES/COEFFS syntax:

% Test when values is passed with periodic boundary conditions:
dom = [0 2*pi];
N = chebop(@(x,u) diff(u,2) + cos(u), dom);
u0 = chebfun('0',dom);
x = chebfun('x',dom);
L = linearize(N, u0, x);
N.bc = 'periodic';
options = cheboppref();
options.discretization = 'values';
isPrefGiven = 1;
out = determineDiscretization(N, L, isPrefGiven, options);
pass(1) = isequal(out.discretization, @trigcolloc);

% Test when coeffs is passed with periodic boundary conditions:
options.discretization = 'coeffs';
isPrefGiven = 1;
out = determineDiscretization(N, L, isPrefGiven, options);
pass(2) = isequal(out.discretization, @ultraS);

% Test when values is passed with periodic boundary conditions and breakpoints:
dom = [0 pi 2*pi];
N = chebop(@(x,u) diff(u,2) + cos(u), dom);
u0 = chebfun('0',dom);
x = chebfun('x',dom);
L = linearize(N, u0, x);
options.discretization = 'values';
isPrefGiven = 1;
out = determineDiscretization(N, L, isPrefGiven, options);
pass(3) = isequal(out.discretization, @chebcolloc2);

% Test when coeffs is passed with periodic boundary conditions and breakpoints:
options.discretization = 'coeffs';
isPrefGiven = 1;
out = determineDiscretization(N, L, isPrefGiven, options);
pass(4) = isequal(out.discretization, @ultraS);

% Test when values is passed with dirichlet boundary conditions:
dom = [-1 1];
N = chebop(@(x,u) diff(u) + exp(u), dom);
u0 = chebfun('0',dom);
x = chebfun('x',dom);
L = linearize(N, u0, x);
N.bc = 'dirichlet';
options.discretization = 'values';
isPrefGiven = 1;
out = determineDiscretization(N, L, isPrefGiven, options);
pass(5) = isequal(out.discretization, @chebcolloc2);

% Test when coeffs is passed with dirichlet boundary conditions:
options.discretization = 'coeffs';
isPrefGiven = 1;
out = determineDiscretization(N, L, isPrefGiven, options);
pass(6) = isequal(out.discretization, @ultraS);

%% Test the default:

% Default with dirichlet:
dom = [-1 1];
N = chebop(@(x,u) diff(u) + sin(u), dom);
u0 = chebfun('0',dom);
x = chebfun('x',dom);
L = linearize(N, u0, x);
N.bc = 'dirichlet';
isPrefGiven = 0; % use the default
out = determineDiscretization(N, L, isPrefGiven, pref); % use pref
pass(7) = isequal(out.discretization, @chebcolloc2);

% Default with periodic:
N.bc = 'periodic';
isPrefGiven = 0; % use the default
out = determineDiscretization(N, L, isPrefGiven, pref); % use pref
pass(8) = isequal(out.discretization, @trigcolloc);

%% Test CHEBCOLLOC1/ULTRAS is passed:

options.discretization = @chebcolloc1;
isPrefGiven = 1; 
out = determineDiscretization(N, L, isPrefGiven, options);
pass(9) = isequal(out.discretization, @chebcolloc1);

options.discretization = @ultraS;
isPrefGiven = 1; 
out = determineDiscretization(N, L, isPrefGiven, options);
pass(10) = isequal(out.discretization, @ultraS);


end