function pass = domainCheck(f, g)
%DOMAINCHECK   True if the domains of two CHEBFUN objects are the same.
%   DOMAINCHECK(F, G) returns TRUE if the endpoints of the domains of the two
%   CHEBFUN objects F and G coincide up to a tolerance depending on their
%   horizontal scales or if both F and G are empty CHEBFUN objects.
%
%   Alternatively, one of F or G may be a two-vector, in which case its values
%   are treated as if they were F.domain([1, end]) or G.domain([1, end]),
%   respectively.
%
% See also HSCALE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The CHEBFUN class uses this function internally to compare the domains of
% CHEBFUN objects before attempting to perform operations on multiple CHEBFUN
% objects that require the CHEBFUN objects to reside on the same interval.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fIsEmpty = isempty(f);
gIsEmpty = isempty(g);

if ( fIsEmpty && gIsEmpty )        % f, g both empty.
    pass = true;
    
elseif ( xor(fIsEmpty, gIsEmpty) ) % Exactly one of f, g is empty.
    pass = false;
    
elseif ( ~isa(g, 'chebfun') )      % f, g both not empty. f is a CHEBFUN.
    hs = hscale(f);
    pass = norm(f(1).domain([1, end]) - g([1, end]), inf) < 1e-15*hs;
    
elseif ( ~isa(f, 'chebfun') )      % f, g both not empty. g is a CHEBFUN.
    hs = hscale(g);
    pass = norm(f([1, end]) - g(1).domain([1, end]), inf) < 1e-15*hs;
    
else                               % f, g both not empty CHEBFUN objects.
    % Grab the hscale:
    hs = max(hscale(f), hscale(g));
    % Compare the domains:
    err = domain(f, 'ends') - domain(g, 'ends');
    % Should be either less than tolerance or NaN (from inf - inf).
    pass = all(abs(err) < 1e-15*hs | isnan(err));
    
end

end
