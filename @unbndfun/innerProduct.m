function out = innerProduct(f, g)
%INNERPRODUCT   Compute the inner product of two UNBNDFUN objects.
%   INNERPRODUCT(F, G) returns the L2 inner product on semi-infinite or doubly-
%   infinite domain of two UNBNDFUN objects. Note that the conjugate of F is
%   taken.
%
%   If F and/or G are array-valued UNBNDFUN objects, then the result is a matrix
%   whose (i,j) entry is the inner product of the ith column of F with the jth
%   column of G.
%
%   The UNBNDFUN objects F and G are assumed to have the same domain. If the
%   domains of F and G are not identical, the output of the method will be
%   garbage though the method doesn't throw a warning or an error message.
%
%   Currently, INNERPRODUCT in UNBNDFUN calls the function with the same name in
%   ONEFUN level. That is,
%
%                        Inf                 1
%                         /                  /
%                        |  ____             |  ____
%   INNERPRODUCT(f, g) = |  f(x) * g(x) dx = |  f(y) * g(y) * m'(y) dy,
%                       /                   /
%                       -Inf                -1
%
%   where m(y) is the forward map of f and m'(y) is the derivative of m(y).
%   Whether the inner product exists depends on the integrability of the
%   integrand conj(f)*g*m', which is handled by the INNERPRODUCT at ONEFUN
%   level.
%
%   Note: UNBNDFUN INNERPRODUCT does not support array-valued inputs.
%
% See also SUM.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Developer note: UNBNDFUN INNERPRODUCT does not support array-valued inputs,
% because SINGFUN does not, and a SINGFUN will be introduced when we multiply by
% the derivative of the nonlinear map.

% Deal with empty case:
if ( isempty(f) || isempty(g) )
    out = [];
    return
end

% If any unrecognizable input, throw an error message.
if ( ~isa(f, 'unbndfun') || ~isa(g, 'unbndfun') )
    error('CHEBFUN:UNBNDFUN:innerProduct:input', ...
        'innerProduct() only operates on two UNBNDFUN objects.');
end

% Assign the output to be the inner product of the onefuns of the input,
% but multiplied by the derivative of the map.

% integrand = conj(f.onefun).*g.onefun.*forDer;
integrand = conj(f).*g;

out = sum(integrand);

end
