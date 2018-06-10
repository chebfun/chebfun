function vout = chebvals2chebvals(vin, kind1, kind2)
%LEGVALS2CHEBVALS   Convert Chebyshev values on first and second-kind grids.
%   V_OUT = CHEBVALS2CHEBVALS(V_IN, 1, 2), converts the vector of V_IN
%   representing values of a polynomial at first-kind Chebyshev points to a
%   vector V_OUT representing the same polynomial evaluated on a second-kind
%   Chebyshev grid.
%
%   V_OUT = CHEBVALS2CHEBVALS(V_IN, 2, 1) is similar except that it maps from a
%   second-kind grid to a first-kind (i.e., p(CHEBPT(N)) --> p(CHEBPT(N,1))).
% 
% See also CHEBVALS2LEGVALS, CHEBPTS, LEGPTS.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

if ( kind1 == kind2 )
    vout = vin;
elseif ( kind1 == 1 && kind2 == 2 )
    c = chebtech1.vals2coeffs(vin);
    vout = chebtech2.coeffs2vals(c);
elseif ( kind1 == 2 && kind2 == 1 )
    c = chebtech2.vals2coeffs(vin);
    vout = chebtech1.coeffs2vals(c);
else
    error('CHEBFUN:chebvals2chebcoeffs:kind', 'Invalid chebkind.');
end

end
