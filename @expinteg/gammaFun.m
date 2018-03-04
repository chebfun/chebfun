function g = gammaFun(j, k)
%GAMMAFUN   Get a function handle to a gamma function.
%   G = GAMMAFUN(J, K) returns a function handle to the gamma function (J, K).
%
% See also EXPINTEG/GAMMAEVAL.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Trivial case:
if ( j == 0 )
    g = @(z) (exp(k*z) - 1)./z;
    
% Compute them recursively using the recurrence formula:
else
    g = @(z) 0*z;
    for m = 1:j
        gm = expinteg.gammaFun(j-m, k);
        g = @(z) g(z) + (-1)^(m-1)/m*gm(z);
    end
    if ( j <= k )
        g = @(z) (g(z) - nchoosek(k, j))./z;
    else
        g = @(z) g(z)./z;
    end
end

end