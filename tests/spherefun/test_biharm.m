function pass = test_biharm( )

tol = 1e4*chebfunpref().cheb2Prefs.chebfun2eps;

%  Test some spherical harmonics
k = 1;
for ell = [1 2 4 5]
    for m = 0:ell
        f = spherefun.sphharm(ell, m);
        lap2 = biharm(f);
        pass(k, 1) = numel(lap2.pivotValues) == numel(f.pivotValues);
        err(k) = SampleError((ell*(ell+1))^2*f, lap2)/(ell*(ell+1))^2;
        pass(k, 2) = SampleError((ell*(ell+1))^2*f, lap2) < (ell*(ell+1))^2*tol;
        k = k+1;
    end
end
pass = pass(:)';

end

function sample_error = SampleError(h, g)
m = 6; 
n = m;
[x, y] = getPoints(m, n);
[L2, T2] = meshgrid(x, y);
F = feval(h, L2, T2);
approx = fevalm(g, x, y);
sample_error = norm(F(:) - approx(:), inf);
end

function [x, y] = getPoints(m, n)

x = trigpts(2*n, [-pi pi]);
y = linspace(0, pi, m).';

end