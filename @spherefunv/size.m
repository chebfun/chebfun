function varargout = size(F, dim) 
%SIZE   Size of a SPHEREFUNV object
%   D = SIZE(F) returns a three-element vector D = [K, M, N]. If F is a column
%   SPHEREFUNV object then K is the number of components in F, N and M are INF.
%   If F is a row vector then K and M are INF and N is the number of components
%   of F.
%
%   [K, M, N] = SIZE(F) returns the dimensions of F as separate output 
%   variables.
%
%   D = SIZE(F, DIM) returns the dimensions specified by the dimension DIM.
%

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(F) ) 
    varargout = {[], []}; 
    return
end

nF = 3; 
tr = F.isTransposed; 

% Manually work out the size: 
if ( ~tr ) 
    K = nF; 
    M = inf; 
    N = inf; 
else
    N = nF;
    K = inf; 
    M = inf; 
end

% Default to vector:
if ( nargin == 1 )
    dim = 0; 
end

% Manually work out what should be displayed:
if ( dim == 1 ) 
    varargout = { K };
elseif ( dim == 2 )
    varargout = { M }; 
elseif ( dim == 3 )
    varargout = { N }; 
elseif ( ( dim == 0 ) && ( nargin == 1 ) )
    if ( nargout <= 1 )
        varargout = { [K, M, N] };
    else
        varargout = { K, M, N }; 
    end
else
    error('SPHEREFUN:SPHEREFUNV:size:dim', 'Unrecognised dimension.');
end

end
