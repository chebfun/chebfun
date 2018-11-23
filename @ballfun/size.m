function varargout = size( F )
%SIZE   Size of a BALLFUN.
%   S = SIZE(F) returns the size of the tensor of expansion coefficients
%   for F, where S = [m,n,p] for an mxnxp tensor of coefficients. 
%   
%   [M, N, P] = SIZE(F) is the same as S = SIZE(F) with S = [M, N, P].

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Grab dimensions of underlying coefficient tensor: 
S = size( F.coeffs );

% If F.coeffs is a matrix, then it has one fiber:  
if ( numel( S ) == 2 )
    S(3) = 1;
end

% Prepare output:
if ( nargout <= 1 ) 
    varargout = { S };
else 
    varargout = { S(1), S(2), S(3) };
end
end
