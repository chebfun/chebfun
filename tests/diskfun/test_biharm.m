function pass = test_biharm( )

tol = 1e7*chebfunpref().cheb2Prefs.chebfun2eps;

%  Test some cylindrical harmonics
%k = [];
%get eigenvalues
err=[];
k=1;
for ell = [ 1 2 3]
    for m = 1:abs(ell)
        %find eigenvalues
        jzero = roots(chebfun(@(x) besselj(ell,x), [sqrt((3/4)^2*pi^2+ell^2) (m+ell/2)*pi]));
        jzero=jzero(m);
        f = diskfun.harmonic(ell, m);
        lap2 = biharm(f);
        %pass(k, 1) = numel(lap2.pivotValues) == numel(f.pivotValues);
        %%NOTE: doesn't always get rank right
        err(k) = SampleError((jzero)^4*f, lap2)/(jzero)^4;
        pass(k) = SampleError((jzero)^4*f, lap2) < 2*(jzero)^4*tol;
        k = k+1;
    end
end


pass = pass(:)';
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



