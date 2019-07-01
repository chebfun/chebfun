function h = times(f, g)
%.*   Pointwise multiplication for BALLFUN.
%   F.*G multiplies BALLFUN F and G. Alternatively F or G could be 
%   a double.
%
% See also MTIMES.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isa(g, 'ballfunv')
    h = times(g,f);
else
    Sf = size(f);
    Sg = size(g);
    S = Sf+Sg;
    F = coeffs3(f,S(1),S(2),S(3));
    G = coeffs3(g,S(1),S(2),S(3));
    Fvals = ballfun.coeffs2vals(F);
    Gvals = ballfun.coeffs2vals(G);
    H = ballfun.vals2coeffs(Fvals.*Gvals);
    h = ballfun(H,'coeffs');
end
end
