function varargout = size( F, varargin )
%SIZE   Size of a BALLFUN.
%   S = SIZE(F) returns the size of the tensor of expansion coefficients
%   for F, where S = [m,n,p] for an mxnxp tensor of coefficients. 
%   
%   [M, N, P] = SIZE(F) is the same as S = SIZE(F) with S = [M, N, P].

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if isempty( F )
    S = [];
else
    % Grab dimensions of underlying coefficient tensor: 
    S = size( F.coeffs );

    % If F.coeffs is a matrix, then it has one fiber:  
    if ( numel( S ) == 2 )
        S(3) = 1;
    end
end

% Prepare output:
if ( nargout <= 1 )
    if nargin > 1
        S = S(varargin{1});
    end
    varargout = { S };
else 
    varargout = { S(1), S(2), S(3) };
end
end
