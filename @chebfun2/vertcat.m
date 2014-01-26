function F = vertcat( varargin )
%VERTCAT Vertical concatenation of chebfun2 objects.
%
% VERTCAT(F, G) is the vertical concatenation of chebfun2 objects F and G, and this
% function returns a chebfun2v object. 
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
            error('CHEBFUN2:VERTCAT', 'Only Chebfun2v objects with 2 or 3 components are valid.');
        else
            Fc = F.components; 
            g = Fc{1};
            h = Fc{2};
            F = chebfun2v( {f, g, h} );
        end
    elseif ( isa(varargin{ 2 }, 'chebfun2' ) )
        % call the chebfun2v constructor.
        F = chebfun2v( varargin );
    end
else
    error('CHEBFUN2:VERTCAT','Cannot vertically concatenate more than three chebfun2 objects.');
end
    

end