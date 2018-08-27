function h = times(f, g)
%.*   BALLFUN multiplication.
%   F.*G multiplies F and G, where F and G may be BALLFUN objects or scalars.
%   If F and/or G is array-valued, the dimensions must match.
%
% See also MTIMES.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Get coeffs: 
F = f.coeffs; 
G = g.coeffs; 

% Pad so F and G are the same size: 
[mf, nf, pf] = size( F ); 
[mg, ng, pg] = size( G );

Fg = zeros(mf+mg, nf+ng, pf+pg); 
Fg(1:mf, 1:nf, 1:pf) = F; 

Gf = zeros(mf+mg, nf+ng, pf+pg); 
Gf(1:mg, 1:ng, 1:pg) = G; 

% Convert to values and do .*: 
Fvals = ballfun.coeffs2vals( Fg );
Gvals = ballfun.coeffs2vals( Gf ); 

% Construct ballfun: 
h = ballfun( Fvals.*Gvals ); 

% Finally simplify ballfun: 
h = simplify( h ); 

end
