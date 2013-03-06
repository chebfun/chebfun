function out = get(f, prop)
%GET    GET method for the CHEBFUN class
%   P = GET(F,PROP) returns the property P specified in the string PROP from
%   the fun F. Valid entries for the string PROP are:
%       'DOMAIN'
%       'FUNS'
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
    case {'lval'}
        out = get(f.funs{1}, 'lval');
    case {'rval'}
        out = get(f.funs{end}, 'rval');
    case {'values', 'coeffs', 'points', 'ishappy', 'hscale', 'vscale', 'epslevel'}
        
        n = numel(f.funs);
        out{n,1} = 0;
        for k = 1:n
            out{k} = get(f.funs{k}, prop);
        end
    otherwise
        error('CHEBFUN:FUNCHEB2:GET:proname', ...
            'Unknown property name ''%s'' for object of type fun.', prop);
end