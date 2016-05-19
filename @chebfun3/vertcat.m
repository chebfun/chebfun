function F = vertcat(varargin)
%VERTCAT   Vertical concatenation of CHEBFUN3 objects.
%   VERTCAT(F, G) is the vertical concatenation of CHEBFUN3 objects F and G. 
%   This function returns a CHEBFUN3V object. 
% 
%   [F; G] is equivalent to VERTCAT(F, G).
%
%   VERTCAT(F) returns the CHEBFUN3 object F. 
% 
% See also CHEBFUN3V.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 ) 
    % VERTCAT of one argument just returns the same thing back to the user:
    F = varargin{1}; 
elseif ( nargin > 1 )
    if ( isa(varargin{2}, 'chebfun3v') )
        f = varargin{1}; 
        F = varargin{2}; 
        if ( F.nComponents > 2 ) 
            error('CHEBFUN:CHEBFUN3:vertcat:tooManyComponents', ...
                'Only CHEBFUN3V objects with 2 or 3 components are valid.');
        else
            Fc = F.components; 
            g = Fc{1};
            h = Fc{2};
            F = chebfun3v({f, g, h});
        end
    elseif ( isa(varargin{ 2 }, 'chebfun3' ) )
        % call the CHEBFUN3V constructor.
        F = chebfun3v( varargin{:} );
    end
else
    error('CHEBFUN:CHEBFUN3:vertcat:tooManyInputs', ...
        'Cannot vertically concatenate more than three CHEBFUN3 objects.');
end
    
end