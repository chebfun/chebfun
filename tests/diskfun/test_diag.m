function pass = test_diag( ) 
% test diskfun diag() command

tol = 10*chebfunpref().cheb2Prefs.chebfun2eps;

% Diagonal of zero function is always zero.
f = diskfun.coeffs2diskfun(0);
pass(1) = iszero(diag(f));

ff = @(th,r) exp(-r.^2.*cos(th).*cos(th));
f = diskfun(ff, 'polar');
rng(10062001);
alp = pi*(1-2*rand(10,1));
for j=1:numel(alp)
    g = chebfun(@(r) ff(alp(j),r), [-1 1]);
    d = diag(f, alp(j));
    pass(1+j) = norm(g-d);
end

try 
    d = diag(f,10*pi);
    pass(2+j) = 0;
catch ME
    pass(2+j) = strcmp(ME.identifier,'CHEBFUN:DISKFUN:diag:angleOutOfRange');
end