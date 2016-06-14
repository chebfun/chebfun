function varargout = length(F)
%LENGTH   Length of a CHEBFUN3.
%   [M, N, P] = LENGTH(F) returns the numbers M, N and P of Chebyshev 
%   (or Fourier) coefficients of the CHEBFUN3 object F in variable X, Y 
%   and Z, respectively. 
%
%   LENGTH(F) returns just the maximum of M, N and P.
%
% See also CHEBFUN3/SIZE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(F) ) 
    % Output:
    if ( nargout <= 1 )
        varargout = {[]};
    else
        varargout = {[], [], []};  
    end
    return
end

if ( iszero(F) ) 
    if ( nargout <= 1 )
        varargout = {1};
    else 
        varargout = {1, 1, 1};
    end
    return
end

m = length(F.cols);
n = length(F.rows);
p = length(F.tubes);

% Output:
if ( nargout <= 1 )
    varargout = {max([m, n, p])};
else
    varargout = {m, n, p};
end

end