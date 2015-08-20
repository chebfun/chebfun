function pass = test_quantumstates(pref)
% Test function for the CHEBFUN method quantumstates().
% AB, 2014/05/08.

%% Setup
dom = [-3 3];   % Domain
x = chebfun(@(x) x, dom);
V = x.^2;       % Potential, harmonic oscillator

%% Default calling sequence
[efuns, evals] = quantumstates(V, 'noplot');

% Extract eigenvalues
D = diag(evals);

% Expect the eigenvalues to go up in steps of .2
err = norm(diff(D) - .2);
pass(1) = ( err < 1e-12);

%% Try a piecewise potential, and a different value for h
h = .24;
V = .1*abs(x-2);
n = 6;

% Solve
[efuns, evals] = quantumstates(V, n, h, 'noplot');

efuns = horzcat(efuns{1:n});

% The Schroedinger operator
op = @(u) -h^2*diff(u,2) + repmat(V, 1, n).*u;

% Did quantumstates() return eigenfunctions and eigenvalues?
err = norm( op(efuns) - efuns*evals );
pass(2) = err < 5e-8;
end
