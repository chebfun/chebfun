function g = fracInt(f, mu, b)
%FRACINT  Fractional integral of a CHEBTECH. 
%   G = FRACINT(F, MU) returns the smooth part, G, of the order MU fractional
%   integral of the CHEBTECH object F. That is, J^{MU}(F(x)) = (1+x)^MU * G.
%   MU is assumed to be in the range (0,1).
%
%   G = FRACINT(F, MU, B) is similar, but returns the smoothpart of the
%   fractional integral of the (1+x)^B*F. That is, J^{MU}((1+x)^B*F(x)) =
%   (1+x)^(MU+B) * G(x). 
%
%   The algorithm used to is to expand F in the Legendre basis and apply the
%   formula [1, (18.17.9)]. The resulting Jacobi polynomial series is converted
%   back to Chebyshev via JAC2CHEB(). In the special case when MU = .5, [1,
%   (18.17.45)] is used instead, as this removes the need for JAC2CHEB().
%
%   References:
%    [1] F.W.J. Olver et al., editors. NIST Handbook of Mathematical Functions.
%    Cambridge University Press, New York, NY, 2010.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
[n, m] = size(f);
if ( m > 1 )
    error('CHEFBFUN:CHEBTECH:fracInt:arrayvalued', ....
        'Array-valued CHEBTECH objects are not supported.');
end
if ( nargin < 3 )
    % Default to b = 0:
    b = 0;
end

if ( b == 0 )
    % Legendre coefficients of f:
    c_leg = cheb2leg(f.coeffs);

    % Determine the new Coefficients via [1, (18.17.9)] or [1, (18.17.45)]:
    if ( mu ~= .5 )
        % Coefficients of the fractional integral:
        c_jac = (c_leg.*beta((0:n-1)'+1, mu))/gamma(mu);
        c_new = jac2cheb(c_jac, -mu, mu);
    else
        % Coefficients of the half-integral:
        z = zeros(1, m);
        scl = (0:n-1)'+.5;
        c_scl = c_leg./(scl*gamma(.5));
        c_new = ( [c_scl ; z] + [z ; c_scl] );

        % For consistence with mu ~= .5, we must divide out by (1+x). (See
        % chebtech/extractBoundaryRoots() for details of this imlpementation.)
        e = ones(n, 1); 
        D = spdiags([.5*e, e, .5*e], 0:2, n, n); 
        D(1) = 1;
        % Compute the new coefficients:
        c_new = [D\c_new(2:end,:) ; z];
    end
    
else
    % P^(0, b) coefficients of f:
    c_jac1 = cheb2jac(f.coeffs, 0, b);

    % Determine the new Coefficients via [1, (18.17.9)]:
    c_jac2 = (c_jac1.*beta((0:n-1)'+b+1, mu))/gamma(mu);
    c_new = jac2cheb(c_jac2, -mu, b+mu);
    
end

% Update f:
g = f;
g.coeffs = c_new;

end
