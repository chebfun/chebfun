function h = times(f, g)
%.*   BALLFUN multiplication.
%   F.*G multiplies F and G, where F and G may be BALLFUN objects or scalars.
%   If F and/or G is array-valued, the dimensions must match.
%
% See also MTIMES.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if (nnz(size(f)-size(g))==0)
    F = f.coeffs;
    G = g.coeffs;
    Fvals = ballfun.coeffs2vals(F);
    Gvals = ballfun.coeffs2vals(G);
    Hvals = Fvals.*Gvals;
    h = ballfun(Hvals);
else
    error('BALLFUN:isequal:unknown', ...
    ['Undefined function ''times'' for different size of ballfun functions : ' ...
     '%s and %s.'], mat2str(size(f)), mat2str(size(g)));
end

end
