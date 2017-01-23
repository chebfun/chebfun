function f = power(f, n)
%.^   CHEBFUN3 power.
%   F.^G returns a CHEBFUN3 F to the scalar power G, a scalar F to the 
%   CHEBFUN3 power G, or a CHEBFUN3 F to the CHEBFUN3 power G. F and/or G 
%   may be complex.
%
%   H = POWER(F, G) is called for the syntax 'F .^ G'.
%
% See also CHEBFUN3/SQRT and CHEBFUN3/COMPOSE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(n) )        % Check for empty objects.
    f = chebfun3();
    
elseif ( isa (f, 'double') && isscalar(f))          % DOUBLE .^ CHEBFUN3    
    op = @(x,y,z) f .^ (feval(n, x, y, z));
    f = chebfun3(op, n.domain);
    
elseif ( isa(n, 'double') && isscalar(n) )          % CHEBFUN3 .^ DOUBLE
   if ( round(n) ~= n )
       % Fractional power --> positive/negative test for the base.
       if ( isreal(f) )
           % Positive/negative test.
           ss = singleSignTest(f);
           if ( ~ss )
               error('CHEBFUN:CHEBFUN3:power:fractional', ...
                   'Sign change detected. Unable to represent the result.'); 
           end
       end
   end
   op = @(x,y,z) feval(f, x, y, z) .^ n;
   f = chebfun3(op, f.domain);

elseif ( (isa(f, 'chebfun3')) && isa(n,'chebfun3') )  % CHEBFUN3 .^ CHEBFUN3
    if ( ~domainCheck(f, n) )          % Check they're on the same domain.
        error('CHEBFUN:CHEBFUN3:power:domain','Domains must be the same.');
    end
    op = @(x,y,z) feval(f, x, y, z) .^ (feval(n, x, y, z));   % Resample
    f = chebfun3(op, f.domain);        % Call constructor
else
    error('CHEBFUN:CHEBFUN3:power:inputs','Inputs must be either CHEBFUN3 objects or scalars.')
end

end