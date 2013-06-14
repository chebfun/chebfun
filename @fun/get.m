function out = get(f, prop)
%GET   GET method for the BNDFUN class.
%   P = GET(F,PROP) returns the property P specified in the string PROP from
%   the fun F. Valid entries for the string PROP are:
%       'DOMAIN' - The domain of F.
%       'MAPPING' - the map used by F.
%       'ONEFUN' - The unmapped representation of F.
%       'VSCALE' - Vertical scale of F.
%       'EPSLEVEL' - Happiness level of F.
%       'LVAL' - Value of F at -1.
%       'RVAL' - Value of F at +1.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

switch prop
    case fieldnames(f)
        % Allow access to any of the properties of F via GET:
        out = f.(prop);
    
    case fieldnames(f.onefun)
        % Allow access to any of the properties of F>onefun via GET:
        out = f.onefun.(prop);

    case {'lval', 'rval'}
        % The value at -1:
        out = get(f.onefun, prop);

    otherwise
        error('CHEBFUN:ONEFUN:GET:proname', ...
            'Unknown property name ''%s'' for object of type ONEFUN.', prop);
end

end