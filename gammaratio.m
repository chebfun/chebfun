function ratio = gammaratio(m, delta)
%GAMMARATIO  Accurately compute a certain ratio of gamma functions.
%   GAMMARATIO(M, D) accurately computes the ratio gamma(M+D)/gamma(M) using a
%   series approximation from [1,p.9]. When M is big and D is small this is more
%   accurate than the naive computation or even expm(gammaln(M+D)-gammaln(M));
%   Reference:
%    [1] N. Hale and A. Townsend, "Fast and accurate computation of Gauss-Legendre 
%        and Gauss-Jacobi quadrature nodes and weights", SIAM J. Sci. Comp., 2013.

if ( m <= 15 || m < delta ) % Determined heuristically
    ratio = expm(gammaln(m+delta)-gammaln(m));
    return
end

% Ensure that 0 < delta < 1:
if ( delta == 0 )
    ratio = 1;
    return
end
fd = floor(delta);
rd = delta - fd;
if ( fd > 1 )
    scl = 1;
    for k = 0:fd-1
        scl = scl*(m+k+rd);
    end
    ratio = scl*gammaratio(m, rd);
    return
end

% Taylor series in expansion 
ds = .5*delta^2/(m-1);
s = ds;
j = 1;
while ( abs(ds/s) > eps/100 && j < 100 ) 
    j = j+1;
    ds = -delta*(j-1)/(j+1)/(m-1)*ds;
    s = s + ds;
end
p2 = exp(s)*sqrt(1+delta/(m-1))*(m-1)^(delta);

% Stirling's series:
g = [1, 1/12, 1/288, -139/51840, -571/2488320, 163879/209018880, ...
    5246819/75246796800, -534703531/902961561600, ...
    -4483131259/86684309913600, 432261921612371/514904800886784000];
f = @(z) sum(g.*[1, cumprod(ones(1, 9)./z)]);

ratio = p2*(f(m+delta-1)/f(m-1));

end
