function pass = domainCheck(f, g)
%DOMAINCHECK   True if the domains of two CHEBFUN objects are the same.
%   DOMAINCHECK(F, G) returns TRUE if the endpoints of the domains of the two
%   CHEBFUNs F and G coincide up to a tolerance depending on their horizontal
%   scales or if both F and G are empty CHEBFUN objects.  
%
% See also HSCALE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The CHEBFUN class uses this function internally to compare the domains of
% CHEBFUN objects before attempting to perform operations on multiple CHEBFUNs
% that require the CHEBFUNs to reside on the same interval.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fIsEmpty = isempty(f);
gIsEmpty = isempty(g);

if ( fIsEmpty && gIsEmpty )        % f, g both empty.
    pass = true;
elseif ( xor(fIsEmpty, gIsEmpty) ) % Exactly one of f, g is empty.
    pass = false;
else                               % f, g both not empty.
    hs = max(hscale(f), hscale(g));    
    pass = norm(f.domain([1, end]) - g.domain([1, end]), inf) < 1e-15*hs;
end

end
