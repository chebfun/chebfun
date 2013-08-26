function out = get(f, prop)
%GET   GET method for the FUN class.
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the FUN object F. Valid entries for the string PROP are:
%       'DOMAIN' - The domain of F.
%       'MAPPING' - The map used by F.
%       'ONEFUN' - The unmapped representation of F on [-1, 1].
%       'VSCALE' - Vertical scale of F.
%       'EPSLEVEL' - Happiness level of F.
%       'LVAL' - Value of F at a (where F.domain = [a,b]).
%       'RVAL' - Value of F at b (where F.domain = [a,b]).

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

switch prop
    case fieldnames(f)
        % Allow access to any of the properties of F via GET:
        out = f.(prop);
    
    case fieldnames(f.onefun)
        % Allow access to any of the properties of F.onefun via GET:
        out = get(f.onefun, prop);

    case {'lval', 'rval'}
        % Get the values at a or b (where f.domain = [a, b]):
        out = get(f.onefun, prop);
        
    case {'points'}
        % Get the underlying grid:
        out = get(f.onefun, prop);

    otherwise
        error('CHEBFUN:ONEFUN:GET:propname', ...
            'Unknown property name ''%s'' for object of type ONEFUN.', prop);
end

end