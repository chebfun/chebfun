function f = power( f, n )
%.^	CHEBFUN2 power. 
%
% F.^G returns a CHEBFUN2 F to the scalar power G, a scalar F to the
% CHEBFUN2 power G, or a CHEBFUN2 F to the CHEBFUN2 power G.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( f ) || isempty( n ) )    % Check for empty objects.
    f = chebfun2();
    
elseif ( isa (f, 'double') )           % double.^CHEBFUN2
    
    op = @(x,y) f .^ ( feval( n, x, y ) );
    f = chebfun2( op, n.domain );
    
elseif ( isa( n, 'double' ) )          % CHEBFUN2.^double

   if ( abs(round( n ) - n) > eps )
       % Positive/negative test.
       [bol, wzero] = singleSignTest(f);
       if ( bol == 0 || wzero == 1 )
           error('CHEBFUN:CHEBFUN2:power:fractional', ...
               'A change of sign/zero has been detected, unable to represent the result.');
       end
   end
    op = @(x,y) feval( f, x, y ) .^ n;
    f = chebfun2( op, f.domain );
    
else                                   % CHEBFUN2.^CHEBFUN2
    
    if ( ~domainCheck(f, n) ) % check they're on the same domain.
        error('CHEBFUN:CHEBFUN2:power:domain','Domains must be the same');
    end
    op = @(x,y) feval( f, x, y ) .^ ( feval( n, x, y ) );
    f = chebfun2( op, f.domain );      % Resample and call constructor.
    
end

end
