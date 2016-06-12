function f = power(f, n)
%.^   CHEBFUN3 power.
%   F.^G returns a CHEBFUN3 F to the scalar power G, a scalar F to the 
%   CHEBFUN3 power G, or a CHEBFUN3 F to the CHEBFUN3 power G. F and/or G 
%   may be complex.
%
%   H = POWER(F, G) is called for the syntax 'F .^ G'.
%
% See also CHEBFUN3/SQRT and CHEBFUN3/COMPOSE.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty(f) || isempty(n) )        % Check for empty objects.
    f = chebfun3();
    
elseif ( isa (f, 'double') )           % double .^ CHEBFUN3    
    if ( ~isscalar(f) )
        error('CHEBFUN:CHEBFUN3:power:scalar','Both inputs must be scalars.')
    else
        op = @(x,y,z) f .^ (feval(n, x, y, z));
        f = chebfun3(op, n.domain);
    end
    
elseif ( isa(n, 'double') )            % CHEBFUN3 .^ double
    if ( ~isscalar(n) )
        error('CHEBFUN:CHEBFUN3:power:scalar','Both inputs must be scalars.')
    elseif ( abs(round(n) - n) > eps )
       % Positive/negative test.
       [bool, wzero] = singleSignTest(f);
       if ( ( bool == 0 ) || ( wzero == 1 ) )
           error('CHEBFUN:CHEBFUN3:power:fractional', ...
               'A change of sign/zero has been detected, unable to represent the result.');
       end
   end
    op = @(x,y,z) feval(f, x, y, z) .^ n;
    f = chebfun3(op, f.domain);
    
else                                   % CHEBFUN3 .^ CHEBFUN3
    if ( ~domainCheck(f, n) )          % Check they're on the same domain.
        error('CHEBFUN:CHEBFUN3:power:domain','Domains must be the same.');
    end
    op = @(x,y,z) feval(f, x, y, z) .^ (feval(n, x, y, z));   % Resample
    f = chebfun3(op, f.domain);        % Call constructor

end

end