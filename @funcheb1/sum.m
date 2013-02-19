function out = sum(f, dim)
%SUM    Definite integral of a FUNCHEB1 on the interval [-1,1].
%   SUM(F) is the integral from -1 to 1 of F.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Get the length of the values:
n = size(f.values, 1);

%%
% Sum across vectorised FUNCHEB1 columns (if dim = 2).
if ( nargin > 1 && dim == 2 )
    f.values = sum(f.values, dim);
    f.coeffs = sum(f.coeffs, dim);
    out = f;
    return
end
   
%%
% Compute integral:

% Trivial cases:
if ( isempty(f) )    % Empty FUNCHEB1
    out = []; 
    return
elseif ( n == 1 )    % Constant FUNCHEB1
    out = 2*f.values;
    return
end

% Integrate in coefficient space:
c = f.coeffs(end:-1:1,:); 
c(2:2:end,:) = 0;

% Note that \int_{-1}^1 T_k(x) dx = \sum_{k = 0}^{n} 2/(1-k^2) for k even, so:
out = [ 2, 0, 2./(1-(2:n-1).^2) ] * c;

