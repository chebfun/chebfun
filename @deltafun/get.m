function out = get(f, prop)
%GET   GET method for the DELTAFUN class.
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the DELTAFUN object F. Valid entries for the string PROP are:
%       'LOCATION'    - ... 
%       'MAGNITUDE'   - ...
%       'FUNPART'     - ...
%       'DOMAIN'      - ...
%       'ISTRANSPOSED - ...

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

switch prop
    case fieldnames(f)
        % Allow access to any of the properties of F via GET:
        out = f.(prop);
    case fieldnames(f.funPart)
        % Access to any of the properties of the smooth part of F:
        out = f.funPart.(prop);
    otherwise
        error('CHEBFUN:DELTAFUN:GET:propname', ...
              'Unknown property name "%s" for object of type DELTAFUN.', prop);
end

end
