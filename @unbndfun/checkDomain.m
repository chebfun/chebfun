function match = checkDomain(f , g)
%CHECKDOMAIN Check whether domains of UNBNDFUN objects F and G match.
%  M = CHECKDOMAIN(F, G) returns 1 if the domains of F and G match, 0 otherwise.

% The tolerance of the match is determined by the hscale of F and G.
tol = max(get(f,'hscale'), get(g,'hscale'))*eps;

% Check whether the difference of the endpoints is less than the tolerance
domf = f.domain;
domg = g.domain;

% If one of f and g has finite domain or both of them have finite domains, then
% return false and bail out.

if all(isfinite(domf)) || all(isfinite(domg))
    match = 0;
    warning('CHEBFUN:UNBNDFUN:checkDomain',...
        'None of the unbndfun objects f and g should have finite domain.');
    return
end

% If both f and g have infinite domains.

if all(isinf(domf)) && all(isinf(domg))
    match = 1;
    return
end

% Now f and g have semi-infinite domains.

if ( isinf(domf(1)) && isinf(domg(1)) ) && ( abs(domf(2) - domg(2)) < tol ) ||...
        ( isinf(domf(2)) && isinf(domg(2)) ) && ( abs(domf(1) - domg(1)) < tol )
    match = 1;
    
else
    match = 0;
    warning('CHEBFUN:UNBNDFUN:checkDomain',...
        'The domain of unbndfun objects f and g do not match.');
end

end