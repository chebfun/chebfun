function out = get(f, prop)
%GET    GET method for the FUN class
%   P = GET(F,PROP) returns the property P specified in the string PROP from
%   the fun F. Valid entries for the string PROP are:
%       'DOMAIN'
%       'MAPPING'
%       'ONEFUN'
%       'VSCALE' - Vertical scale of F.
%       'ISHAPPY' - Is F happy?
%       'EPSLEVEL' - Happiness level of F.
%       'POINTS' - Grid corresponding to F.
%       'LVAL'
%       'RVAL'

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

switch prop
    case fieldnames(f)
        % Allow access to any of F's properties via GET.
        out = f.(prop);
    case 'points'
         out = f.mapping.for(f.onefun.(prop));
    case {'lval', 'rval', 'values', 'coeffs', 'ishappy', 'hscale', 'vscale', 'epslevel'}
        out = get(f.onefun, prop);        
    otherwise
        error('CHEBFUN:FUN:GET:proname', ...
            'Unknown property name ''%s'' for object of type fun.', prop);
end