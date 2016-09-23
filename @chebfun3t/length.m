function varargout = length(f)
%LENGTH   Length of a CHEBFUN3T object.
%   This returns size of tensor of Chebyshev coefficients. If just one 
%   output has been asked, then it is the max number of coefficients in 
%   three dimensions.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) ) 
    % Output:
    if ( nargout <= 1 )
        varargout = {[]};
    else
        varargout = {[], [], []};  
    end
    return
end

coeffs = f.coeffs;

if ( ndims(coeffs) < 3 )                % If the function is univariate or
    [m,n] = size(coeffs);               % bivariate, artificially put 1 as
    p = 1;                              % the size in the 3rd dimension.
else
    [m,n,p] = size(coeffs);
end
    
% Output:
if ( nargout <= 1 )
    varargout = {max([m, n, p])};
else
    varargout = {m, n, p};  
end

end