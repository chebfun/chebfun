function f = power( f, n )
%.^	Chebfun2 power. 
%
% F.^G returns a chebfun2 F to the scalar power G, a scalar F to the
% chebfun2 power G, or a chebfun2 F to the chebfun2 power G.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty( f ) || isempty( n ) )    % check for empty objects.
    f = chebfun2;
    return;
end

if ( isa (f, 'double') )                % double.^chebfun2
    
    op = @(x,y) f .^ ( feval( n, x, y ) );
    f = chebfun2( op, n.domain );
    
elseif ( isa( n, 'double' ) )            % chebfun2.^double

   if ( abs(round( n ) - n) > eps )
       % positive/negative test.
       [bol, wzero] = chebfun2.singleSignTest(f);
       if ( bol == 0 || wzero == 1 )
           error('CHEBFUN2:POWER:FRACTIONAL',...
               'A change of sign/zero has been detected, unable to represent the result.');
       end
   end
    op = @(x,y) feval( f, x, y ) .^ n;
    f = chebfun2( op, f.domain );
else                                  % chebfun2.^chebfun2
    if ( ~chebfun2.domainCheck(f, n) ) % check they're on the same domain.
        error('CHEBFUN2:power:domain','Domains must be the same');
    end
    
    op = @(x,y) feval( f, x, y ) .^ ( feval( n, x, y ) );
    f = chebfun2( op, f.domain );       % resample and call constructor.
end

end