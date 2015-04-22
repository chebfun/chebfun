function f = power( f, n )
%.^	     Pointwise power of a LOWRANKAPPROX. 
%
% F.^G returns a LOWRANKAPPROX F to the scalar power G, a scalar F to the
% CHEBFUN2 power G, or a LOWRANKAPPROX F to the LOWRANKAPPROX power G.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isempty( f ) || isempty( n ) )    % Check for empty objects.
    
    f = chebfun2( );
    
elseif ( isa (f, 'double') )           % double.^LOWRANKAPPROX
    
    op = @(x,y) f .^ ( feval( n, x, y ) );
    f = chebfun2( op, n.domain );
    
elseif ( isa( n, 'double' ) )          % LOWRANKAPPROX.^double

   if ( abs(round( n ) - n) > eps )
       % Positive/negative test.
       [bol, wzero] = singleSignTest(f);
       if ( bol == 0 || wzero == 1 )
           error('CHEBFUN:LOWRANKAPPROX:power:fractional', ...
               'A change of sign/zero has been detected, unable to represent the result.');
       end
   end
    op = @(x,y) feval( f, x, y ) .^ n;
    f = chebfun2( op, f.domain );
    
else                                   % LOWRANKAPPROX.^LOWRANKAPPROX
    
    if ( ~domainCheck(f, n) ) % check they're on the same domain.
        error('CHEBFUN:LOWRANKAPPROX:power:domain','Domains must be the same');
    end
    op = @(x,y) feval( f, x, y ) .^ ( feval( n, x, y ) );
    f = chebfun2( op, f.domain );      % Resample and call constructor.
    
end

end
