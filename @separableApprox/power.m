function f = power(f, n)
%.^   Pointwise power of a SEPARABLEAPPROX.
%
%   F.^G returns a SEPARABLEAPPROX F to the scalar power G, a scalar F to the
%   SEPARABLEAPPROX power G, or a SEPARABLEAPPROX F to the SEPARABLEAPPROX
%   power G.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) )    % Check for empty objects.

    return
    
elseif ( isempty(n) )
    
    f = n; 
    return
    
elseif ( isa(f, 'double') )           % double.^SEPARABLEAPPROX
    
    f = compose(f, @power, n);
    
elseif ( isa( n, 'double' ) )          % SEPARABLEAPPROX.^double

   if ( abs(round(n) - n) > eps )
       if ( isreal(f) )
           % Positive/negative test.
           bool = singleSignTest(f);  % Returns TRUE if there is no sign change.
           if ( ~bool )
               error('CHEBFUN:SEPARABLEAPPROX:power:fractional', ...
                   'Sign change detected. Unable to represent the result.'); 
           end
       end
   end
   
   % Call compose of child object:
   f = compose(f, @power , n);
    
else                                   % SEPARABLEAPPROX.^SEPARABLEAPPROX
    
    if ( ~domainCheck(f, n) ) % Check they're on the same domain.
        error('CHEBFUN:SEPARABLEAPPROX:power:domain','Domains must be the same');
    end
    f = compose(f, @power, n);      % Resample and call constructor.
    
end

end