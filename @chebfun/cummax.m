function g = cummax(f)
%CUMMAX   Cumulative maximum of a CHEBFUN.
%   G = CUMMAX(F) is the cumulative maximum of a row or column CHEBFUN F
%   over its domain of definition.
%
% Example:
%
%   rng(0), f = randnfun(.01); plot(f,'b',cummax(f),'r')
%
% See also CUMMIN.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Trivial case:
if ( isempty(f) )
    return
end

[rows, cols] = size(f);
% Quasimatrices are not supported:
if ( (rows == Inf && cols > 1) || (cols == Inf && rows > 1) )
    error('CHEBFUN:CHEBFUN:cummax:quasi', ...
        'CUMMAX does not currently support quasimatrices.');
end

% Check if we are dealing with a row CHEBFUN
rowChebfun = 0;
if rows == 1
    rowChebfun = 1;
    f = f';
end
% Extract the local maxima of F
[fx, x] = max(f,'local');
% Perform a discrete CUMMAX on the local maxima of F 
xCum = x(1);
fxCum = fx(1);
for i = 2:length(fx)
    if ( fx(i) > fxCum(end) )
        fxCum = [fxCum, fx(i)];
        xCum = [xCum, x(i)];
    end
end
x = xCum;
fx = fxCum;

% construct the CUMMAX of F
if( x(1) > f.domain(1) )
    fi = simplify(restrict(f, [f.domain(1), x(1)]));
    funs = fi.funs;
else
    funs = {};
end

for i = 1:length(x)-1
    fi = max(simplify(restrict(f, [x(i), x(i+1)])), fx(i));
    funs = [funs, fi.funs];
end

if( x(end) < f.domain(end) )
    fi = chebfun(fx(end), [x(end), f.domain(end)]);
    funs = [funs, fi.funs];
end

% Assemble the new CHEBFUN:
g = chebfun(funs);
if rowChebfun
    g = g';
end

end