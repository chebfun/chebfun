function out = get(f, prop)
%GET   GET method for the SINGFUN class.
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the SINGFUN object F. Valid entries for the string PROP are:
%       'EXPONENTS'  - The exponents of F at the end points -1 and 1.
%       'SMOOTHPART' - The smooth part of F on [-1, 1], which is a CHEBTECH.
%   GET also allows access to the properties of the smoothPart.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

switch prop
    
    case fieldnames(f.smoothPart)
        % Access to any of the properties of the smooth part of F:
        out = f.smoothPart.(prop);
        
    case fieldnames(f)
        % Allow access to any of the properties of F via GET:
        out = f.(prop);
        
    otherwise
        error('CHEBFUN:SINGFUN:GET:propname', ...
            'Unknown property name "%s" for object of type SINGFUN.', prop);
        
end

end
