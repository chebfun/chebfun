function h = mtimes(f, g)
%*   BALLFUN multiplication.
%   a*F and F*a multiplies the BALLFUN F by the scalar a.
%
% See also TIMES.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

fIsBallfun = isa(f, 'ballfun');
gIsBallfun = isa(g, 'ballfun');
fIsBallfunv = isa(g, 'ballfunv');

if (fIsBallfun && isnumeric(g))
    h = ballfun(g*f.coeffs, 'coeffs');
elseif (isnumeric(f) && gIsBallfun)
    h = ballfun(f*g.coeffs, 'coeffs');
elseif (fIsBallfun && gIsBallfun)
    h = times(f,g);
elseif fIsBallfunv
    h = mtimes(g, f);
else
    error('BALLFUN:mtimes:unknown', ...
          ['Undefined function ''mtimes'' for input arguments of type ' ...
           '%s and %s.'], class(f), class(g));
end
end
