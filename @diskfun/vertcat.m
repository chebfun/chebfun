function F = vertcat(varargin)
%VERTCAT   Vertical concatenation of DISKFUN objects.
%   K = VERTCAT(F, G, H) is the vertical concatenation of DISKFUN objects F, 
%   G, and H. The output K is a DISKFUNV object.
% 
%   [F ; G] is equivalent to VERTCAT(F, G).
%
%   VERTCAT(F) returns the DISKFUN F. 
% 
% See also DISKFUNV.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


if ( nargin == 1 ) 
    % VERTCAT of one argument just returns the same thing back to the user:
    F = varargin{1};
    
elseif ( nargin == 2 )
    if ( all(cellfun(@(F) isa(F,'diskfun'), varargin)) )
        F = diskfunv(varargin{:});
        
    else
        error('DISKFUN:vertcat:tooManyComponents', ...
            'Only DISKFUN objects are valid to concatenate.');
    end
    
else
    error('DISKFUN:vertcat:tooManyInputs', ...
        'Can only vertically concatenate two DISKFUN objects.');
end
    
end