function F = vertcat( varargin )
%VERTCAT Vertical concatenation of CHEBFUN2 objects.
%
% VERTCAT(F, G) is the vertical concatenation of CHEBFUN2 objects F and G, and this
% function returns a CHEBFUN2V object. 
% 
% [F ; G] is different syntax for VERTCAT(F, G)
% 
% See also CHEBFUN2V.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( nargin > 1 )
    if ( isa( varargin{ 2 } , 'chebfun2v') )
        f = varargin{1}; 
        F = varargin{2}; 
        if ( F.nComponents > 2 ) 
            error('CHEBFUN2:VERTCAT', 'Only CHEBFUN2V objects with 2 or 3 components are valid.');
        else
            Fc = F.components; 
            g = Fc{1};
            h = Fc{2};
            F = chebfun2v( {f, g, h} );
        end
    elseif ( isa(varargin{ 2 }, 'chebfun2' ) )
        % call the CHEBFUN2V constructor.
        F = chebfun2v( varargin );
    end
else
    error('CHEBFUN2:VERTCAT','Cannot vertically concatenate more than three CHEBFUN2 objects.');
end
    

end