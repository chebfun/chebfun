function g = min( f, varargin )
%MIN  Minimum value of a chebfun2
% 
% MIN(F), returns a chebfun representing the minimum of the chebfun2 along
% the y direction, i.e 
% 
%   MIN(F) = @(x) min( F ( x, : ) ) 
% 
% MIN(F,[],DIM) returns a chebfun representing the minimum of f along the
% DIM direction.
% 
% See also MAX. 

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check: 
if ( isempty( f ) ) 
    error('CHEBFUN2:MIN:INPUT','Chebfun2 is empty');
end

% Default to maximum along the y direction:
if ( nargin == 1 )   
    dim=1;
end

if ( nargin > 1 )   % do not allow min(F,G). 
    if ( ~isempty(varargin{1}) )
        error('CHEBFUN2:MIN','Unable to minimise two chebfun2 objects.');
    end
end

if ( nargin == 2 )  % default to maximum along the y direction. 
    dim = 1;
end

if ( nargin == 3 )  % what direction do we want to max in? 
    dim = varargin{ 2 };
end


% min(f) = - max ( -f )  
g = -max( -f, [], dim );

end