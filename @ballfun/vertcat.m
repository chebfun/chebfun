function F = vertcat( varargin ) 
%VERTCAT   Vertical concatenation of BALLFUN objects.
%
% K = VERTCAT(F, G, H) is the vertical concatenation of BALLFUN objects F, 
% G, and H. The function K is a BALLFUNV object. 
% 
% [F ; G ; H] is equivalent to VERTCAT(F, G, H).
%
% VERTCAT(F, G) returns an error. BALLFUNV objects have three components.
%
% VERTCAT(F) returns the BALLFUN F. 
% 
% See also BALLFUNV.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 ) 
    % VERTCAT of one argument just returns the same thing back to the user:
    F = varargin{1}; 
elseif ( nargin == 3 )
    if all( cellfun(@(F) isa(F,'ballfun'), varargin) )
        F = ballfunv( varargin{:} );
    else
        error('BALLFUN:vertcat:tooManyComponents', ...
            'Only BALLFUN objects are valid to concatenate.');
    end
else
    error('BALLFUN:vertcat:tooManyInputs', ...
        'Can only vertically concatenate three BALLFUN objects.');
end
    
end