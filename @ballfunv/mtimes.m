function h = mtimes(f, g)
%*   BALLFUNV multiplication.
%   A*F and F*A multiplies the BALLFUNV F by the scalar A.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

fIsBallfunv = isa(f, 'ballfunv');
gIsBallfunv = isa(g, 'ballfunv');

if (fIsBallfunv && isnumeric(g))
    F = f.comp;
    h = ballfunv(g*F{1},g*F{2},g*F{3});
elseif (isnumeric(f) && gIsBallfunv)
    G = g.comp;
    h = ballfunv(f*G{1},f*G{2},f*G{3});
else 
    error('BALLFUNV:mtimes:unknown', ...
          ['Undefined function ''mtimes'' for input arguments of type ' ...
           '%s and %s.'], class(f), class(g));
end
end
