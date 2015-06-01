function [f, g] = overlap(f, g)
%OVERLAP   Overlap the domain of two CHEBFUN objects.
%   [FOUT, GOUT] = OVERLAP(F, G) returns two CHEBFUN objects FOUT and GOUT such
%   that DOMAIN(FOUT) == DOMAIN(GOUT) and F(x) = FOUT(x), G(x) = GOUT(x) for all
%   x in the domain of F and G. If F and/or G are have more than one column/row
%   then all columns of FOUT and GOUT will have the same domain.
%
% See also RESTRICT.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check that the domains are valid:
if ( ~domainCheck(f, g) )
    error('CHEBFUN:CHEBFUN:overlap:domains', ...
        'Inconsistent domains; domain(f) ~= domain(g).')
end

% Grab the domains:
fDom = domain(f);
gDom = domain(g);

if ( (numel(f) == 1) && (numel(g) == 1) )
    if ( (numel(fDom) == numel(gDom)) && all(fDom == gDom) )
        % Trivial case: Nothing to do!
        return
    end
    [f, g] = tweakDomain(f, g);
    fDom = domain(f);
    gDom = domain(g);
else
    g = tweakDomain(g, fDom);
    gDom = domain(g);
    f = tweakDomain(f, gDom);
    fDom = domain(f);
end

% Obtain the domain for the output:
newDom = union(fDom, gDom);

% Restrict f and g:
f = restrict(f, newDom);
g = restrict(g, newDom);

end
