function g = power(f,n)
%.^   BALLFUN power.
%   F.^G returns a BALLFUN F to the scalar power G.
%
%   H = POWER(F, G) is called for the syntax 'F .^ G'.
%
% See also SQRT.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Fast power implementation
if (n >= 0)
    S = size(f);
    g = ballfun(@(r,lam,th)1,S);
    while (n>0)
        % If power is odd
        if (mod(n,2) == 1)
            g = g*f;
        end
        n = floor(n/2);
        f = f*f;
    end 
else
    error('BALLFUN:isequal:unknown', ...
    ['Undefined function ''power'' for negative power : ' ...
     '%d.'], n);
end
end
