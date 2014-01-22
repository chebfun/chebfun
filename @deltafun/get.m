function out = get(f, prop)
%GET   GET method for the DELTAFUN class.
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the DELTAFUN object F. Valid entries for the string PROP are:
%
%       'LOCATION'    - Location of the delta functions 
%       'DELTAMAG'    - Magnitude of the delta functions
%       'FUNPART'     - The smooth function contained in DELTAFUN.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

switch prop
    case fieldnames(f)
        % Allow access to any of the properties of F via GET:
        out = f.(prop);
    case fieldnames(f.funPart)
        % Access to any of the properties of the smooth part of F:
        out = f.funPart.(prop);
    case fieldnames(f.funPart.onefun)
        out = f.funPart.onefun.(prop);        
    case {'lval', 'rval'}
        % Get the values at a or b (where f.domain = [a, b]):
        out = get(f.funPart, prop);
        
    case {'points'}
        % Get the underlying grid:
        out = get(f.funPart, prop);

    otherwise
        error('DELTAFUN:GET:propname', ...
              'Unknown property name "%s" for object of type DELTAFUN.', prop);
end

end
