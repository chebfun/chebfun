function F = vertcat(varargin)
%VERTCAT   Vertical concatenation of SPHEREFUN objects.
%
% VERTCAT(F, G) is the vertical concatenation of SPHEREFUN objects F and G.
% This function returns a SPHEREFUNV object. 
% 
% [F; G] is equivalent to VERTCAT(F, G).
%
% VERTCAT(F) returns the SPHEREFUN F. 
% 
% See also SPHEREFUNV.

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