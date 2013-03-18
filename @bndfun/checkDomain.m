function match = checkDomain(f , g)
%CHECKDOMAIN Check whether domains of BNDFUN objects F and G match.
%   M = CHECKDOMAIN(F, G) returns 1 if the domains of F and G match, 0 otherwise.

% TODO: Do we really need a checkDomain method at the bndfun level?

% The tolerance of the match is determined by the hscale of F and G.
tol = max(get(f,'hscale'), get(g,'hscale'))*eps;

% Check whether the difference of the endpoints is less than the tolerance
match = all(abs(f.domain - g.domain) < tol);

end