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

m = mf + mg;
n = nf + ng;
p = pf + pg;

Fg = zeros(m,n,p); 
Fg(1:mf,floor(n/2)+1-floor(nf/2):floor(n/2)+nf-floor(nf/2),floor(p/2)+1-floor(pf/2):floor(p/2)+pf-floor(pf/2)) = F; 

Gf = zeros(m,n,p); 
Gf(1:mg,floor(n/2)+1-floor(ng/2):floor(n/2)+ng-floor(ng/2),floor(p/2)+1-floor(pg/2):floor(p/2)+pg-floor(pg/2)) = G; 

% Convert to values and do .*: 
Fvals = ballfun.coeffs2vals( Fg );
Gvals = ballfun.coeffs2vals( Gf ); 

% Construct ballfun: 
h = ballfun( Fvals.*Gvals ); 

% Finally simplify ballfun: 
h = simplify( h ); 

end
