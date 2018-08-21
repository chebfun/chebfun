function h = mtimes(f, g)
%*   BALLFUN multiplication.
%   A*F and F*A multiplies the BALLFUN F by the scalar A.
%
% See also TIMES.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

fIsBallfun = isa(f, 'ballfun');
gIsBallfun = isa(g, 'ballfun');

if (fIsBallfun && isnumeric(g))
    H = g*f.coeffs;
    h = ballfun(H);
elseif (isnumeric(f) && gIsBallfun)
    H = f*g.coeffs;
    h = ballfun(H);
elseif (fIsBallfun && gIsBallfun)
    h = times(f,g);
else
    error('BALLFUN:mtimes:unknown', ...
          ['Undefined function ''mtimes'' for input arguments of type ' ...
           '%s and %s.'], class(f), class(g));
end

end
