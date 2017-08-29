function varargout = size( F, dim ) 
%SIZE   Size of a DISKFUNV object
%   D = SIZE(F) returns a three-element vector D = [K,M,N]. If F is a column
%   DISKFUNV object then K is the number of components in F, N and M are INF.
%   If F is a row vector then K and M are INF and N is the number of 
%   components of F.
%
%   [K,M,N] = SIZE(F) returns the dimensions of F as separate output 
%   variables.
%
%   D = SIZE(F,DIM) returns the dimensions specified by the dimension DIM.
%
% See also CHEBFUN2/SIZE. 

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information. 

% Empty check: 
if ( isempty( F ) ) 
    varargout = { [], [] }; 
    return
end

% Get number of components: 
nF = F.nComponents; 
trans = F.isTransposed; 

% Manually work out the size: 
if ( ~trans ) 
    K = nF; 
    M = inf; 
    N = inf; 
else
    N = nF;
    K = inf; 
    M = inf; 
end

% Default to a vector:
if ( nargin == 1 )
    dim = 0; 
end

% Manually work out what should be displayed:
if ( dim == 1 ) 
    % K = SIZE( F, 1 ): 
    varargout = { K };
    
elseif ( dim == 2 )
    % M = SIZE( F, 2 ):
    varargout = { M }; 
    
elseif ( dim == 3 )
    % N = SIZE( F, 3 ):
    varargout = { N }; 
    
elseif ( ( dim == 0 ) && ( nargin == 1 ) )
    
    if ( nargout <= 1 )
        % SIZE(F) and K = SIZE( F ): 
        varargout = { [K, M, N] };
        
    else
        % [K, M, N] = SIZE( F ): 
        varargout = { K, M, N }; 
        
    end
    
else
    error('CHEBFUN:DISKFUNV:size:dim', 'Unrecognised dimension.');
end

end