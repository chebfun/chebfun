function pass = test_laplacian( )


tol = 4*1e4*chebfunpref().cheb2Prefs.chebfun2eps;

%  Test some cylindrical harmonics
%k = [];
%get eigenvalues
err=[];
k=1;
for ell = [ 1 2 3 4]
    for m = 1:abs(ell)
        %find eigenvalues
        jzero = roots(chebfun(@(x) besselj(ell,x), [sqrt((3/4)^2*pi^2+ell^2) (m+ell/2)*pi]));
        jzero=jzero(m);
        f = diskfun.harmonic(ell, m);
        lapf = laplacian(f);
        %pass(k, 1) = numel(lap2.pivotValues) == numel(f.pivotValues);
        %%NOTE: doesn't always get rank right
        err(k) = SampleError(-(jzero)^2*f, lapf)/(jzero)^2;
        pass(k) = SampleError(-(jzero)^2*f, lapf) < (jzero)^2*tol;
        k = k+1;
    end
end

%exponential
g = @(x,y) exp(-3*x).*y.^4; 
f = diskfun(g);
lapf = laplacian(f); 
exact = @(x,y) (-3)^2*exp(-3*x).*y.^4+12*y.^2.*exp(-3*x);
pass(11) = SampleError(diskfun(exact), lapf) < 4e1*tol;

%check that lap gives the same result
lapf = lap(f); 
pass(12) = SampleError(diskfun(exact), lapf) < 4e1*tol;

end

function sample_error = SampleError(h, g)
m = 7; 
n = m-1;
[x, y] = getPoints(m, n);
[L2, T2] = meshgrid(x, y);
F = feval(h, L2, T2, 'polar');
approx = fevalm(g, x, y);
sample_error = norm(F(:) - approx(:), inf);
end

function [x, y] = getPoints(m, n)

x = trigpts(2*n, [-pi pi]);
y = chebpts(m);
y = y(ceil(m/2):end); 

end






