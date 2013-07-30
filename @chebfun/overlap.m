function [f, g] = overlap(f, g)
%OVERLAP    Overlap the domain of two chebfun objects.
%   [FOUT, GOUT] = OVERLAP(F ,G) returns two chebfuns such that FOUT.domain ==
%   GOUT.domain and F(x) = FOUT(x), G(x) = GOUT(x) for all x in domain of F.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check that the domains are valid:
if ( ~domainCheck(f, g) )
    error('CHEBFUN:overlap:domains', ...
        'Inconsistent domains; domain(f) ~= domain(g).')
end

% Obtain the domains of f and g:
fDom = f.domain;
gDom = g.domain;

% Take the union of the two domains: (Neither fDom or gDom will be empty.)
newDom = union(fDom, gDom);

% If f and g are both empty, there is nothing to do.
if ( isempty(newDom) )
    return
end

if ( length(fDom) ~= length(gDom) || ~all(fDom == gDom) )
    % Breakpoints do not match. We have work to do.
    % [TODO]: Should we allow a tolerance so as not to introcude tiny intervals?

    % Compute the new objects using RESTRICT():
    f = restrict(f, newDom);
    g = restrict(g, newDom);

end

% Sort out the impulses: (Pad to have the same length)
fImps = f.impulses;
gImps = g.impulses;
fRows = size(fImps, 1);
gRows = size(gImps, 1);
maxRows = max(fRows, gRows);
if ( maxRows ~= fRows )
    f.impulses = [ fImps ; zeros(maxRows - fRows, length(newDom))];
end
if ( maxRows ~= gRows )
    g.impulses = [ gImps ; zeros(maxRows - gRows, length(newDom))];
end

end