function f = set(f, prop, propName)
%SET     Set a DISKFUNV property.
%
%   SET(F, PROP, PROPNAME) sets the property PROP of the DISKFUNV F to PROP.
%   Valid options for PROP are
%
%   'coords' - The coordinate system used to evaluate F. This can be set to
%   either 'polar' (for polar coordinates) or 'cart' (for Cartesian
%   coordinates).
%
% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( F ) )
    return
end

switch prop
    
    case 'coords'
        
        % Set DISKFUNV coordinate system on each component:
        f.coords = propName;
        f.components{1}.coords = propName;
        f.components{2}.coords = propName;
        
    otherwise
        error('CHEBFUN:DISKFUN:set:Prop', 'Unknown DISKFUN property.')
        
end
end