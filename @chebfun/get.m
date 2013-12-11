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

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% [TODO]: Include a get(f, 'numCols') ( = size(f.funs{1}, 2) if f is not empty).

switch prop
    case {'domain', 'ends'}
        out = domain(f);
        
    case fieldnames(f)
        % Allow access to any of F's properties via GET.
        if ( numel(f) > 1 )
            % TODO: Is this the right thing to do?
            warning('CHEBFUN:get:quasi', ...
                'CHEBFUN/GET() returns only the property information of the first column.')
        end
        out = f(1).(prop);
        
    case 'lval'
        out = cell(1, numColumns(f));
        for k = 1:numel(f)
            out{k} = get(f(k).funs{1}, 'lval');
        end
        out = cell2mat(out);
        if ( f(1).isTransposed )
            out = out.';
        end
            
    case 'rval'
        out = cell(1, numColumns(f));
        for k = 1:numel(f)
            out{k} = get(f(k).funs{end}, 'rval');
        end
        out = cell2mat(out);
        if ( f(1).isTransposed )
            out = out.';
        end
        
    case {'lval-local', 'rval-local'}
        if ( iszero(f) )
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
        for j = 1:numel(f)
            for k = 1:numel(f(j).funs)
                out = all(get(f(j).funs{k}, 'ishappy'));
                if ( ~out ), return, end
            end
        end
        
    case 'hscale'
        out = hscale(f);   

    case 'vscale'
        out = vscale(f); 

    case 'epslevel'
        out = epslevel(f);

    case {'values', 'coeffs', 'points', 'epslevel-local', 'vscale-local', 'hscale-local'}
        % Adjust the PROP variable in the case of -local:
        idx = strfind(prop, '-');
        if ( ~isempty(idx) )
            prop = prop(1:idx-1);
        end
        
        % Initialise the output:
        m = numel(f);
        n = 0;
        for j = 1:m
            n = max(n, numel(f(j).funs));
        end
        out = cell(n, m);
        
        % Loop over each of the columns and each of the funs:
        for k = 1:m
            for j = 1:numel(f(k).funs)
                out{j,k} = get(f(k).funs{j}, prop);
            end
        end
        if ( ~isempty(idx) && min(m, n) == 1 )
            out = cell2mat(out);
        end
        if ( f(1).isTransposed )
            out = out.';
        end
        
    otherwise
        error('CHEBFUN:CHEBFUN:GET:propname', ...
            'Unknown property name ''%s'' for object of type CHEBFUN.', prop);
        
end
