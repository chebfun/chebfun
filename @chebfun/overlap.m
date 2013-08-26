function [f, g] = overlap(f, g)
%OVERLAP   Overlap the domain of two CHEBFUN objects.
%   [FOUT, GOUT] = OVERLAP(F ,G) returns two CHEBFUNs such that FOUT.DOMAIN ==
%   GOUT.DOMAIN and F(x) = FOUT(x), G(x) = GOUT(x) for all x in the  domain of
%   F. The third dimension of the IMPULSES fields of F and G will be padded with
%   zeros if necessary so that FOUT.IMPULSES and GOUT.IMPULSES store impulse
%   data to the same order.

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

% Pad the impulse arrays so that outputs store impulse data to the same order.
fImps = f.impulses;
gImps = g.impulses;
maxImpOrder = max(size(fImps, 3), size(gImps, 3));
padF = zeros(size(fImps, 1), size(fImps, 2), maxImpOrder - size(fImps, 3));
f.impulses = cat(3, fImps, padF);
padG = zeros(size(gImps, 1), size(gImps, 2), maxImpOrder - size(gImps, 3));
g.impulses = cat(3, gImps, padG);

end
