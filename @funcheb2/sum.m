function out = sum(f, dim)
%SUM    Definite integral of a FUNCHEB2 on the interval [-1,1].
%   SUM(F) is the integral from -1 to 1 of F.
%
%   If F is a vector-valued FUNCHEB2 then the result is a row vector containing
%   the definite integrals of each column.
%
%   SUM(F, 2) sums over the second dimension of F, i.e., adds each of its
%   columns. If F is a scalar-valued FUNCHEB2, this simply returns F.
%
% See also CUMSUM, DIFF.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Get the length of the values:
n = size(f.values, 1);

%%
% Sum across vectorised FUNCHEB2 columns if dim = 2:
if ( nargin > 1 && dim == 2 )
    f.values = sum(f.values, dim);
    f.coeffs = sum(f.coeffs, dim);
    out = f;
    return
end
   
%%
% Compute the integral if dim = 1:

% Trivial cases:
if ( isempty(f) )    % Empty FUNCHEB2
    out = []; 
    return
elseif ( n == 1 )    % Constant FUNCHEB2
    out = 2*f.values;
    return
end

% Evaluate the integral by using the Chebyshev coefficients (see Thm. 19.2 of
% Trefethen, Approximation Theory and Approximation Practice, SIAM, 2013, which
% states that \int_{-1}^1 T_k(x) dx = 2/(1-k^2) for k even):
c = f.coeffs(end:-1:1,:); 
c(2:2:end,:) = 0;
out = [ 2, 0, 2./(1-(2:n-1).^2) ] * c;

end
