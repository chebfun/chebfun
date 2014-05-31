function [f, g] = overlap(f, g)
%OVERLAP   Overlap the domain of two CHEBFUN objects.
%   [FOUT, GOUT] = OVERLAP(F, G) returns two CHEBFUNs such that FOUT.DOMAIN ==
%   GOUT.DOMAIN and F(x) = FOUT(x), G(x) = GOUT(x) for all x in the domain of
%   F. 
%
%   If both F and/or G are have more than one column/row then either both must
%   have the same number of columns/rows or one must have only a single column,
%   otherwise an error is thrown. In the case where, say, G has a single column
%   and F has many, the breakpoints in all columns of F will be unified.
%
% See also RESTRICT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check that the domains are valid:
if ( ~domainCheck(f, g) )
    error('CHEBFUN:overlap:domains', ...
        'Inconsistent domains; domain(f) ~= domain(g).')
end
% Check the columns are valid:
if ( numColumns(f) ~= numColumns(g) &&  ...
    numColumns(f) ~= 1 && numColumns(g) ~= 1 )
    error('CHEBFUN:overlap:dim', 'Matrix dimensions must agree.');
end

if ( numel(f) == 1 && numel(g) == 1)
    % CHEBFUN - CHEBFUN
    [f, g] = columnOverlap(f, g);
elseif ( numel(f) > 1 && numel(g) > 1)
    % QUASIMATRIX - QUASIMATRIX
    for k = 1:numel(f)
        [f(k), g(k)] = columnOverlap(f(k), g(k));
    end
    
elseif ( numel(f) > 1 )
    % QUASIMATRIX - CHEBFUN
    
    if ( numColumns(g) == 1)
        % Scalar expansion in g:
        [f(1), g] = columnOverlap(f(1), g); % Make breaks in f(1) the same as g.
        f = restrict(f, domain(f));         % Make breaks in f the same as f(1).
        [f(1), g] = columnOverlap(f(1), g); % Make breaks in g the same as in f.
    else
        gCell = mat2cell(g);
        g = 0*f;
        for k = 1:numel(f)
            [f(k), g(k)] = columnOverlap(f(k), gCell{k});
        end
    end
    
else % if ( numel(g) > 1 )    
    % CHEBFUN - QUASIMATRIX

    if ( numColumns(f) == 1)
        % Scalar expansion in f (see above):
        [f, g(1)] = columnOverlap(f, g(1));
        g = restrict(g, domain(g));
        [f, g(1)] = columnOverlap(f, g(1));
    else
        fCell = mat2cell(f);
        f = 0*g;
        for k = 1:numel(f)
            [f(k), g(k)] = columnOverlap(fCell{k}, g(k));
        end
    end
end    

end

function [f, g] = columnOverlap(f, g)

% If f and g are both empty, there is nothing to do:
if ( isempty(f) && isempty(g) )
    return
end

if ( (length(f.domain) ~= length(g.domain)) || ~all(f.domain == g.domain) )
    % Tweak the domain to prevent the introduction of tiny intervals:
    [f, g] = tweakDomain(f, g);

    % Take the union of the two domains: (NB: At least one of fDom or gDom is
    % nonempty, so we don't need to worry about the orientation of the output of
    % union().)
    newDomain = union(f.domain, g.domain);
    
    % Breakpoints do not match. Compute the new objects using RESTRICT():
    f = restrict(f, newDomain);
    g = restrict(g, newDomain);
    
end

end
