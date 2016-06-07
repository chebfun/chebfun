function F = vertcat( varargin )
%VERTCAT   Vertical concatenation of SPHEREFUN objects.
%
% K = VERTCAT(F, G, H) is the vertical concatenation of SPHEREFUN objects F, 
% G, and H. The function K is a SPHEREFUNV object. 
% 
% [F ; G ; H] is equivalent to VERTCAT(F, G, H).
%
% VERTCAT(F, G) returns an error. SPHEREFUNV objects have three components.
%
% VERTCAT(F) returns the SPHEREFUN F. 
% 
% See also SPHEREFUNV.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


if ( nargin == 1 ) 
    % VERTCAT of one argument just returns the same thing back to the user:
    F = varargin{1}; 
elseif ( nargin == 3 )
    if all( cellfun(@(F) isa(F,'spherefun'), varargin) )
        F = spherefunv( varargin{:} );
    else
        error('SPHEREFUN:vertcat:tooManyComponents', ...
            'Only SPHEREFUN objects are valid to concatenate.');
    end
else
    error('SPHEREFUN:vertcat:tooManyInputs', ...
        'Can only vertically concatenate three SPHEREFUN objects.');
end
    
end