function f = power(f, n)
%.^   CHEBFUN3T power.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(n) )    % Check for empty objects.
    f = chebfun3t();
    
elseif ( isa (f, 'double') )        % double .^ CHEBFUN3T
    op = @(x,y,z) f .^ (feval(n, x, y, z));
    f = chebfun3t(op, n.domain);
    
elseif ( isa(n, 'double') )          % CHEBFUN3T .^ double
    op = @(x,y,z) feval(f, x, y, z) .^ n;
    f = chebfun3t(op, f.domain);
    
else                                  % CHEBFUN3T .^ CHEBFUN3T
    
    if ( ~domainCheck(f, n) ) % check they're on the same domain.
        error('CHEBFUN:CHEBFUN3T:power:domain','Domains must be the same');
    end
    op = @(x,y,z) feval(f, x, y, z) .^ (feval( n, x, y, z));
    f = chebfun3t(op, f.domain);     % Resample and call constructor.
end

end