function out = get(f, prop)
%GET    GET method for the ONEFUN class
%   P = GET(F,PROP) returns the property P specified in the string PROP from
%   the fun F. Valid entries for the string PROP are:
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
    case {'points', 'lval', 'rval', 'ishappy', 'vscale', 'epslevel'}
        out = get(f.smoothfun, prop);
    otherwise
        error('CHEBFUN:ONEFUN:GET:proname', ...
            'Unknown property name ''%s'' for object of type onefun.', prop);
end