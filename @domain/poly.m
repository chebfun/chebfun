function f = poly(v, d)
%POLY   Convert roots to CHEBFUN.
%   POLY(V, D), when V is a columns vector and D is a domain, returns a CHEBFUN
%   of degree length(V) and domain D whose roots are the elements of V.
%
% See also CHEBFUN/POLY.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Deal with array input:
if ( min(size(v)) > 1 )
    f = cell(1, size(v, 2));
    for k = 1:size(v, 2)
        f{k} = poly(v(:,k), d);
    end
    f = horzcat(f{:});
    return
end

d = double(d);
N = length(v);

% Remove infs, and return NaN if NaN present:
v = v(~isinf(v));
if ( any(isnan(v)) )
    f = chebfun(NaN, d);
    return
end

% Return empty CHEBFUN if v is empty:
if ( N == 0 )
    f = chebfun();
    return
end

% Leja ordering:
[ignored, j] = max(abs(v));
z = v(j);
v(j) = [];
for k = 1:(N-1)
    P = zeros(N - k, 1);
    for l = 1:(N-k)
        P(l) = prod(z - v(l));
    end
    [ignored, j] = max(abs(P));
    z(k+1) = v(j);
    v(j) = [];
end
v = z;

% Evaluate at Chebyshev points:
x = chebpts(N+1, d);
p = ones(N+1, 1);
for k = 1:N
    p = p.*(x - v(k));
end

% Contruct the CHEBFUN:
f = chebfun(p, d, 'tech', @chebtech2);

end
