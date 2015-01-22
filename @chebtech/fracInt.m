function g = fracInt(f, mu)
%FRACINT  Fractional integral of a CHEBTECH. 
%   G = FRACINT(F, MU) returns the smooth part, G, of the order MU fractional
%   integral of the CHEBTECH object F. That is, (J^{MU}F)(x) = (1+x)^MU * G.
%   MU is assumed to be in the range (0,1).
%
%   The algorithm used to is to expand F in the Legendre basis and apply the
%   formula [1, (18.17.9)]. The resulting Jacobi polynomial series is converted
%   back to Chebyshev via JAC2CHEB(). In the special case when MU = .5, [1,
%   (18.17.45)] is used instead, as this removes the need for JAC2CHEB().
%
%   References:
%    [1] F.W.J. Olver et al., editors. NIST Handbook of Mathematical Functions.
%    Cambridge University Press, New York, NY, 2010.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
[n, m] = size(f);
if ( m > 1 )
    error('CHEFBFUN:CHEBTECH:fracInt:arrayvalued', ....
        'Array-valued CHEBTECH objects are not supported.');
end
if ( ~isa(f, 'chebtech') )
    error('CHEFBFUN:CHEBTECH:halfInt:chebtech', ....
        'FRACINT() only supports CHEBTECH objects.');
end

% Legendre coefficients of f:
c_leg = cheb2leg(f.coeffs);

% Determine the new coefricients via [1, (18.17.9)] or [1, (18.17.45)]:
if ( mu ~= .5 )
    % Coefficients of the fractional integral:
    c_jac = (c_leg.*beta((0:n-1)'+1, mu))/gamma(mu);
    c_new = jac2cheb(c_jac, -mu, mu);
else
    % Coefficients of the half-integral:
    z = zeros(1, m);
    scl = (0:n-1)'+.5;
    c_scl = c_leg./scl/gamma(.5);
    c_new = ( [c_scl ; z] + [z ; c_scl] );

    % For consistence with mu ~= .5, we must divide out by (1+x). (See
    % chebtech/extractBoundaryRoots() for details of this imlpementation.)
    e = ones(n, 1); 
    D = spdiags([.5*e, e, .5*e], 0:2, n, n); 
    D(1) = 1;
    % Compute the new coefficients:
    c_new = [D\c_new(2:end,:) ; z];
end

% Update f:
g = f;
g.coeffs = c_new;
g.vscale = getvscl(f);
g.epslevel = updateEpslevel(f);

end