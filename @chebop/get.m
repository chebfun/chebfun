function out = get(N, prop)
%GET   GET method for the CHEBTECH2 class.
%   P = GET(N, PROP) returns the property P specified in the string PROP from
%   the CHEBOP N.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

switch prop
    case fieldnames(N)
        % Allow access to any of the properties of F via GET:
        out = N.(prop);
    otherwise
        error('CHEBFUN:CHEBOP:GET:proname', ...
            'Unknown property name ''%s'' for object of type CHEBOP.', prop);
end

end
