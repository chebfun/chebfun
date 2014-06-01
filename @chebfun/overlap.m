function [f, g] = overlap(f, g)
%OVERLAP   Overlap the domain of two CHEBFUN objects.
%   [FOUT, GOUT] = OVERLAP(F, G) returns two CHEBFUN objects FOUT and GOUT such
%   that DOMAIN(FOUT) == DOMAIN(GOUT) and F(x) = FOUT(x), G(x) = GOUT(x) for all
%   x in the domain of F and G. If F and/or G are have more than one column/row
%   then all columns of FOUT and GOUT will have the same domain.
%
% See also RESTRICT.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check that the domains are valid:
if ( ~domainCheck(f, g) )
    error('CHEBFUN:overlap:domains', ...
        'Inconsistent domains; domain(f) ~= domain(g).')
end

% Obtain the domain for the output:
fDom = domain(f);
g = tweakDomain(g, fDom);
gDom = domain(g);
f = tweakDomain(f, gDom);
fDom = domain(f);
newDom = union(fDom, gDom);

% Restrict f and g:
f = restrict(f, newDom);
g = restrict(g, newDom);

end
