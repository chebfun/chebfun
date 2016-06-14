function F = vertcat(varargin)
%VERTCAT   Vertical concatenation of CHEBFUN3 objects.
%   VERTCAT(G, H, W) is the vertical concatenation of CHEBFUN3 objects F, G
%   and H. This function returns a CHEBFUN3V object.
% 
%   [G; H; W] is equivalent to VERTCAT(G, H, W).
%
%   VERTCAT(G) returns the CHEBFUN3 object G.
% 
% See also CHEBFUN3V.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1 ) 
    % VERTCAT of one argument just returns the same thing back to the user:
    F = varargin{1};
elseif ( nargin > 1 && nargin <= 3)
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
        
    elseif ( isa(varargin{2}, 'chebfun3' ) )
        % Call the CHEBFUN3V constructor:
        F = chebfun3v(varargin{:});
    else
        error('CHEBFUN:CHEBFUN3:vertcat:badInputs', ...
            'Inputs must be either CHEBFUN3 or CHEBFUN3V objects.');
    end
    
else
    error('CHEBFUN:CHEBFUN3:vertcat:tooManyInputs', ...
        'Cannot vertically concatenate more than three objects.');
end

end