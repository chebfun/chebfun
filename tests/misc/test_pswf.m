% Test file for pswf and pswfpts.m.

function pass = test_pswf(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

%% Test PSWF:
P = pswf(0:9, 4);

% Check orthogonality:
pass(1) = norm(P'*P - diag(2./(1:2:19))) < 1e-10;

% Check number of roots:
pass(2) = numel(roots(P(:,end))) == 9;

% Check a precomputed value:
P4pi05 = -0.335179227182412994123178747; % From Wolfram Alpha
pass(3) = abs(feval(pswf(4,pi), 0.5) - P4pi05) < 1e-10;

%% Test PSWFPTS:

% Check nodes are roots of P10:
[x, w] = pswfpts(9, 4);
pass(4) = norm(x - sort(roots(P(:,end)))) < 1e-10;

% Check integral of P(1:10):
pass(5) = norm(sum(P) - w*P(x,:)) < 1e-10;

% Check integral of cosine(k*x):
k = 51;
f = @(x) cos(k*x);
I = 2*sin(k)/k;
[x, w] = pswfpts(k, k);
err = abs(w*f(x) - I);
pass(6) = err < 1e-10;

%% Test GGQ:

% Check nodes are roots of P10:
[x, w] = pswfpts(5, 4, [-1 1], 'ggq');

% CHeck integral of P(1:10):
pass(7) = norm(sum(P) - w*P(x,:)) < 1e-10;

% Check integral of cosine(k*x):
k = 51;
f = @(x) cos(k*x);
I = 2*sin(k)/k;
[x, w] = pswfpts(k, k, [-1 1], 'ggq');
err = abs(w*f(x) - I);
pass(8) = err < 1e-10;

end