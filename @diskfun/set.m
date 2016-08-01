function f = set(f, prop, propName)
%SET: currently can only be used to change setting of the 'coords' 
% parameter.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.
 
    switch prop
        case 'coords'
                f.coords = propName;
        otherwise
            error('CHEBFUN:DISKFUN:set:Prop', ...
                'Unknown DISKFUN property.')
    end
    
    

end
