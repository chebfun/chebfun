function pass = test_bcsyntax

% Test whether the various syntaxes for BC settings are valid. They are NOT
% checked for correctness.

% TAD, 12 May 2014

warnState = warning('off', 'CHEBFUN:CHEBOP:parseBC:keywordbc');

N = chebop(@(u) diff(u,2) - exp(u));
N.lbc = 1;
N.rbc = @(u) diff(u)-2;

N = chebop(@(u) diff(u,2) - exp(u));
N.bc = @(u) [ u(0)-1, sum(u) ];

N = chebop(@(u) diff(u,3) - exp(u));
N.lbc = 'dirichlet';
N.rbc = {1,'neumann'};
N.bc = @(u) u(0)-1;

pass = true;

warning(warnState)

end

