function h = mtimes(f, g)
%*  mtimes for BALLFUNV.
%
%  c*F or F*c multiplies each component of the BALLFUNV F by the scalar c.
%
%   See also TIMES.

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( f ) || isempty( g )
    h = ballfunv();
    return
end

fIsBallfunv = isa(f, 'ballfunv');
gIsBallfunv = isa(g, 'ballfunv');
fIsBallfun = isa(f, 'ballfun');
gIsBallfun = isa(g, 'ballfun');

if ((isnumeric(g) || gIsBallfun) && fIsBallfunv)
    F = f.comp;
    h = ballfunv(g*F{1},g*F{2},g*F{3});
elseif ((isnumeric(f) || fIsBallfun) && gIsBallfunv)
    G = g.comp;
    h = ballfunv(f*G{1},f*G{2},f*G{3});
elseif fIsBallfunv && gIsBallfunv
    h = f.*g;
else 
    error('BALLFUNV:mtimes:unknown', ...
          ['Undefined function ''mtimes'' for input arguments of type ' ...
           '%s and %s.'], class(f), class(g));
end
end
