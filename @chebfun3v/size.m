function varargout = size(F, dim)
%SIZE   size of a CHEBFUN3V object.
%   D = SIZE(F) returns a four-element vector D = [K, M, N, P]. If F is a 
%   column CHEBFUN3V object then K is the number of components in F, and 
%   M, N and P are INF. If F is a row vector then K, M and N are INF and P 
%   is the number of components in F.
%
%   [K, M, N, P] = SIZE(F) returns the dimensions of F in separate outputs.
%
%   D = SIZE(F, DIM) returns the dimension specified by DIM.
%
% See also CHEBFUN3/SIZE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( F ) ) 
    varargout = { [], [], [], [] };
    return
end

nF = F.nComponents; 
tr = F.isTransposed; 

% Manually work out the size: 
if ( ~tr ) 
    K = nF; 
    M = Inf; 
    N = Inf; 
    P = Inf;
else
    P = nF;
    K = Inf; 
    M = Inf;
    N = Inf;
end

% Default to vector:
if ( nargin == 1 )
    dim = 0; 
end

% Manually work out what should be displayed:
if ( dim == 1 ) 
    varargout = {K};
elseif ( dim == 2 )
    varargout = {M}; 
elseif ( dim == 3 )
    varargout = {N}; 
elseif ( ( dim == 0 ) && ( nargin == 1 ) )
    if ( nargout <= 1 )
        varargout = {[K, M, N, P]};
    else
        varargout = {K, M, N, P};
    end
else
    error('CHEBFUN:CHEBFUN3V:size:dim', 'Unrecognised dimension.');
end

end