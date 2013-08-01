function out = get(f, prop)
%GET    GET method for the CHEBFUN class
%   P = GET(F,PROP) returns the property P specified in the string PROP from
%   the fun F. Valid entries for the string PROP are:
%       'DOMAIN'         - The domain of definintion of F.
%       'FUNS'           - The piecewise smooth components of F.
%       'VSCALE'         - Vertical scale of F.
%       'VSCALE-LOCAL'   - Local vertical scales of F.
%       'HSCALE'         - Horizontal scale of F.
%       'HSCALE-LOCAL'   - Local horizontal scales of F.
%       'ISHAPPY'        - Is F happy?
%       'EPSLEVEL'       - Approximate accuracy estimate of F.
%       'EPSLEVEL-LOCAL' - Approximate accuracy estimate of F's components.
%       'POINTS'         - Grid corresponding to F.
%       'LVAL'           - Value(s) of F at lefthand side of domain.
%       'RVAL'           - Value(s) of F at righthand side of domain.

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
    case 'ishappy'
        n = numel(f.funs);
        out(n,1) = 0;
        for k = 1:n
            out(k) = all(get(f.funs{k}, prop));
        end
        out = all(out);
    case 'hscale'
        out = hscale(f);                
    case 'hscale-local'
        n = numel(f.funs);
        out(n,1) = 0;
        for k = 1:n
            out(k) = get(f.funs{k}, 'hscale');
        end   
    case 'vscale'
        out = vscale(f);        
    case 'vscale-local'
        n = numel(f.funs);
        m = size(f.funs{1}, 2);
        out = zeros(n, m);
        for k = 1:n
            out(k,:) = get(f.funs{k}, 'vscale');
        end
    case 'epslevel'
        out = epslevel(f);
    case 'epslevel-local'
        n = numel(f.funs);
        out(n,1) = 0;
        for k = 1:n
            out(k) = get(f.funs{k}, 'epslevel');
        end            
    case {'values', 'coeffs', 'points'}
        n = numel(f.funs);
        out{n,1} = 0;
        for k = 1:n
            out{k} = get(f.funs{k}, prop);
        end
    case {'ends'}
        out = f.domain;
    otherwise
        error('CHEBFUN:CHEBFUN:GET:proname', ...
            'Unknown property name ''%s'' for object of type chebfun.', prop);
end
