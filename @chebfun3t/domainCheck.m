function out = domainCheck(f, g)
%DOMAINCHECK   True if the domains of two CHEBFUN3T objects are the same.
%   DOMAINCHECK(F, G) returns TRUE if the domains of the two CHEBFUN3T 
%   objects F and G coincide up to a tolerance depending on their
%   horizontal scales or if both F and G are empty CHEBFUN objects.
%
% See also CHEBFUN/DOMAINCHECK.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The CHEBFUN3T class uses this function internally to compare the domains of
% CHEBFUN3T objects before attempting to perform operations on multiple
% CHEBFUN3T objects that require the CHEBFUN3T objects to reside on the same interval.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Empty check: 
if ( isempty(f) && isempty(g) )
    out = true;
    return
elseif ( xor(isempty(f), isempty(g) ) )
    out = false;
    return
end

% Compute INF norm of the domains:
hscaleF = norm(f.domain, inf);
hscaleG = norm(g.domain, inf);
hscale = max(hscaleF, hscaleG);

% Compare the domains:
err = f.domain - g.domain;

% Should be less than tolerance.
out = all(abs(err) < 1e-15*hscale);

end