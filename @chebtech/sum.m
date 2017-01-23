function out = sum(f, dim)
%SUM   Definite integral of a CHEBTECH on the interval [-1,1].
%   SUM(F) is the integral of F from -1 to 1.
%
%   If F is an array-valued CHEBTECH, then the result is a row vector
%   containing the definite integrals of each column.
%
%   SUM(F, 2) sums over the second dimension of F, i.e., adds up its columns.
%   If F is a scalar-valued CHEBTECH, this simply returns F.
%
% See also CUMSUM, DIFF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% Get the length of the values:
n = size(f.coeffs, 1);

%%
% Sum across array-valued CHEBTECH columns if dim = 2:
if ( nargin > 1 && dim == 2 )
    f.coeffs = sum(f.coeffs, dim);
    out = f;
    return
end
   
%%
% Compute the integral if dim = 1:

% Trivial cases:
if ( isempty(f) )    % Empty CHEBTECH
    out = []; 
    return
elseif ( n == 1 )    % Constant CHEBTECH
    out = 2*f.coeffs;
    return
end

% Evaluate the integral by using the Chebyshev coefficients (see Thm. 19.2 of
% Trefethen, Approximation Theory and Approximation Practice, SIAM, 2013, which
% states that \int_{-1}^1 T_k(x) dx = 2/(1-k^2) for k even):
c = f.coeffs; 
c(2:2:end,:) = 0;
out = [ 2, 0, 2./(1-(2:n-1).^2) ] * c;

end
