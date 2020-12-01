% Test file for pswfpts.m.

function pass = test_pswfpts(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

P = pswf(0:9, 4);

% Check nodes are roots of P9:
[x, w] = pswfpts(9, 4);
pass(1) = norm(x - sort(roots(P(:,end)))) < 1e-10;

% Check integrals of P(0:9):
pass(2) = norm(sum(P) - w*P(x,:)) < 1e-10;

% Check integral of cosine(k*x):
k = 51;
f = @(x) cos(k*x);
I = 2*sin(k)/k;
[x, w] = pswfpts(k, k);
err = abs(w*f(x) - I);
pass(3) = err < 1e-10;

%% Test GGQ:

% Check nodes are roots of P9:
[x, w] = pswfpts(5, 4, [-1 1], 'ggq');

% Check integral of P(0:9):
pass(4) = norm(sum(P) - w*P(x,:)) < 1e-10;

% Check integral of cosine(k*x):
k = 51;
f = @(x) cos(k*x);
I = 2*sin(k)/k;
[x, w] = pswfpts(k, k, [-1 1], 'ggq');
err = abs(w*f(x) - I);
pass(5) = err < 1e-10;

% Check symmetry:
[x,w] = pswfpts(10,pi,[-7,7]);
pass(6) = (norm(w-fliplr(w)) == 0) & (norm(x+flipud(x)) == 0);
[x,w] = pswfpts(15,sqrt(2),[-7,7],'ggq');
pass(7) = (norm(w-fliplr(w)) == 0) & (norm(x+flipud(x)) == 0);

% Check c > N regime:
[x, w] = pswfpts(21, 80, 'ggq');
P = pswf(40,80);
pass(8) = norm(sum(P) - w*P(x,:)) < 1e-10;

end
