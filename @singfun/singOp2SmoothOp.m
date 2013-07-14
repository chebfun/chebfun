function op = singOp2SmoothOp( op, exponents, tol )
%SINGOP2SMOOTHOP   Converts the original functiona handle OP, which is
%   presumably singular at the end points to smooth operator by factoring
%   out the singular factors
%
% See also SINGFUN

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( all(exponents < -100*tol) )
    % both exponents are non trivial
    op = @(x) op(x)./((1+x).^(exponents(1)).*(1-x).^(exponents(2)));
elseif ( exponents(1) < -100*tol )
    % singular factor at the left end point only
    op = @(x) op(x)./(1+x).^(exponents(1));
elseif ( exponents(2) < -100*tol )
    % singular factor at the right end point only
    op = @(x) op(x)./(1-x).^(exponents(2));
end