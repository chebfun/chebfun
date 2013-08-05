function [f, g] = overlap(f, g)
%OVERLAP   Overlap the domain of two CHEBFUN objects.
%   [FOUT, GOUT] = OVERLAP(F ,G) returns two CHEBFUNs such that FOUT.DOMAIN ==
%   GOUT.DOMAIN and F(x) = FOUT(x), G(x) = GOUT(x) for all x in domain of F.
%   Additionally, the third dimension of the IMPULSES fields of F and G will
%   be padded with zeros if necessary so that FOUT.IMPULSES and GOUT.IMPULSES
%   store impulse data to the same order.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check that the domains are valid:
if ( ~domainCheck(f, g) )
    error('CHEBFUN:overlap:domains', ...
        'Inconsistent domains; domain(f) ~= domain(g).')
end

% If f and g are both empty, there is nothing to do:
if ( isempty(f) && isempty(g) )
    return
end

% Obtain the domains of f and g:
fDom = f.domain;
gDom = g.domain;

% Take the union of the two domains: (NB:  At least one of fDom or gDom is
% nonempty, so we don't need to worry about the orientation of the output of
% union().)
newDom = union(fDom, gDom);

if ( (length(fDom) ~= length(gDom)) || ~all(fDom == gDom) )
    % Breakpoints do not match. We have work to do.
    % [TODO]: Should we allow a tolerance so as not to introduce tiny intervals?

    % Compute the new objects using RESTRICT():
    f = restrict(f, newDom);
    g = restrict(g, newDom);

end

% Pad the impulse arrays so that outputs store impulse data to the same order.
fImps = f.impulses;
gImps = g.impulses;
fMaxImpOrder = size(fImps, 3);
gMaxImpOrder = size(gImps, 3);
maxImpOrder = max(fMaxImpOrder, gMaxImpOrder);
if ( maxImpOrder ~= fMaxImpOrder )
    nFuns = length(newDom);
    nCols = size(f, 2);
    f.impulses = cat(3, fImps, zeros(nFuns, nCols, maxImpOrder - fMaxImpOrder));
end
if ( maxImpOrder ~= gMaxImpOrder )
    nFuns = length(newDom);
    nCols = size(g, 2);
    g.impulses = cat(3, gImps, zeros(nFuns, nCols, maxImpOrder - gMaxImpOrder));
end

end
