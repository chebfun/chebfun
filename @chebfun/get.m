function out = get(f, prop)
%GET   GET method for the CHEBFUN class
%   P = GET(F, PROP) returns the property P specified in the string PROP from
%   the CHEBFUN F. Valid entries for the string PROP are:
%       'DOMAIN'         - The domain of definintion of F.
%       'FUNS'           - The piecewise smooth components of F.
%       'VSCALE'         - Vertical scale of F.
%       'VSCALE-LOCAL'   - Local vertical scales of F.
%       'HSCALE'         - Horizontal scale of F.
%       'HSCALE-LOCAL'   - Local horizontal scales of F.
%       'ISHAPPY'        - Is F happy?
%       'EPSLEVEL'       - Approximate accuracy estimate of F.
%       'EPSLEVEL-LOCAL' - Approximate accuracy estimate of F's components.
%       'LVAL'           - Value(s) of F at lefthand side of domain.
%       'RVAL'           - Value(s) of F at righthand side of domain.
%       'LVAL-LOCAL      - Value(s) of F's FUNs at left sides of their domains.
%       'RVAL-LOCAL'     - Value(s) of F's FUNs at right sides of their domains.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: Include a get(f, 'numCols') ( = size(f.funs{1}, 2) if f is not empty).

switch prop
    case fieldnames(f)
        % Allow access to any of F's properties via GET.
        out = f.(prop);
    case 'lval'
        out = get(f.funs{1}, 'lval');
        if ( f.isTransposed )
            out = out.';
        end
    case 'rval'
        out = get(f.funs{end}, 'rval');
        if ( f.isTransposed )
            out = out.';
        end
    case {'lval-local', 'rval-local'}
        if ( isempty(f) )
            out = [];
            return
        end
        lrval = prop(1:4);
        numFuns = numel(f.funs);
        numCols = size(f.funs{1}, 2);
        out = zeros(numFuns, numCols);
        for k = 1:numFuns
            out(k,:) = get(f.funs{k}, lrval);
        end
        if ( f.isTransposed )
            out = out.';
        end        
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
        m = min(size(f));
        n = numel(f.funs);
        out = zeros(n, m);
        for k = 1:n
            out(k,:) = get(f.funs{k}, 'vscale');
        end
    case 'epslevel'
        out = epslevel(f);
    case 'epslevel-local'
        n = numel(f.funs);
        numCols = min(size(f));
        % TODO:
%         numCols = numColumns(f);
        out(n,numCols) = 0;
        for k = 1:n
            out(k,:) = get(f.funs{k}, 'epslevel');
        end            
    case {'values', 'coeffs', 'points'}
        n = numel(f.funs);
        out = cell(n, 1);
        for k = 1:n
            out{k} = get(f.funs{k}, prop);
        end
        if (numel(out) == 1)
            out = out{:};
        end
    case 'ends'
        out = f.domain;
    otherwise
        error('CHEBFUN:CHEBFUN:GET:propname', ...
            'Unknown property name ''%s'' for object of type CHEBFUN.', prop);
end
