function c = ultra2ultra(c, lam_in, lam_out)
%ULTRA2ULTRA   Convert between Ultraspherical (US) expansions.
%   C_OUT = ULTRA2ULTRA(C_IN, LAM_IN, LAM_OUT) converts the vector C_IN of
%   US-LAM_IN coefficients to a vector of US-LAM_OUT coefficients.
%
%   This code is essentially a wrapper for the JAC2JAC code based on [1].
%
%   References:
%     [1] A. Townsend, M. Webb, and S. Olver, "Fast polynomial transforms
%     based on Toeplitz and Hankel matrices", submitted, 2016.
%
% See also JAC2JAC.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

n = length(c) - 1;

% Scale input from US to Jacobi basis:
c = c./scl(lam_in, n);
% Call JAC2JAC:
c = jac2jac(c, lam_in-.5, lam_in-.5, lam_out-.5, lam_out-.5);
% Scale output from Jacobi to US:
c = c.*scl(lam_out, n);

end

function s = scl(lam, n)
% Scaling from Jacobi to US polynomials. See DLMF Table 18.3.1.
    if ( lam == 0 )
        nn = (0:n-1).';
        s = [1 ; cumprod((nn+.5)./(nn+1))];
    else
        nn = (0:n).';
        s = ( gamma(2*lam) ./ gamma(lam+.5) ) * ...
            exp( gammaln(lam+.5+nn) - gammaln(2*lam+nn) );
    end
end