function f = power(f, n)
%.^   CHEBFUN3 power.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(n) )    % Check for empty objects.
    f = chebfun3();
    
elseif ( isa (f, 'double') )           % double.^CHEBFUN3    
    op = @(x,y,z) f .^ (feval(n, x, y, z));
    f = chebfun3(op, n.domain, 'fiberDim', 3);
    
elseif ( isa(n, 'double') )          % CHEBFUN3.^double
   if ( abs(round(n) - n) > eps )
       % Positive/negative test.
       [bool, wzero] = singleSignTest(f);
       if ( bool == 0 || wzero == 1 )
           error('CHEBFUN:CHEBFUN3:power:fractional', ...
               'A change of sign/zero has been detected, unable to represent the result.');
       end
   end
    op = @(x,y,z) feval(f, x, y, z) .^ n;
    f = chebfun3(op, f.domain, 'fiberDim', 3);
    
else                                   % CHEBFUN3.^CHEBFUN3
    if ( ~domainCheck(f, n) ) % check they're on the same domain.
        error('CHEBFUN:CHEBFUN3:power:domain','Domains must be the same');
    end
    op = @(x,y,z) feval(f, x, y, z) .^ (feval(n, x, y, z));
    f = chebfun3(op, f.domain, 'fiberDim', 3);  % Resample and call constructor.

end

end