function [f, g] = overlap(f, g)
%OVERLAP   Overlap the domain of two CHEBFUN objects.
%   [FOUT, GOUT] = OVERLAP(F, G) returns two CHEBFUNs such that FOUT.DOMAIN ==
%   GOUT.DOMAIN and F(x) = FOUT(x), G(x) = GOUT(x) for all x in the domain of
%   F. 
%
%   If F and G are array-valued, they must have the same numer of columns/rows,
%   else an error is thrown.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Should we allow scalar expansion?

% Check that the domains are valid:
if ( ~domainCheck(f, g) )
    error('CHEBFUN:overlap:domains', ...
        'Inconsistent domains; domain(f) ~= domain(g).')
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
    gCell = mat2cell(g);
    if ( numel(gCell) ~= numel(f) )
        error('CHEBFUN:overlap:dim', 'Matrix dimensions must agree.')
    else
        g = 0*f;
        for k = 1:numel(f)
            [f(k), g(k)] = columnOverlap(f(k), gCell{k});
        end
    end
       
else % if ( numel(g) > 1 )    
    % CHEBFUN - QUASIMATRIX
    fCell = mat2cell(f);
    if ( numel(fCell) ~= numel(g) )
        error('CHEBFUN:overlap:dim', 'Matrix dimensions must agree.')
    else
        f = 0*g;
        for k = 1:numel(g)
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
