function h = rdivide(f,g)
%./   Pointwise right divide of SEPARABLEAPPROX objects.
%   F./G if F is a SEPARABLEAPPROX and G is a double this returns (1/G)*F
%
%   F./G if F is a double and G is a v this returns F/G, but this does
%   not work if G becomes numerically close to zero.
%
%   F./G we do not allow F and G to both be SEPARABLEAPPROX object.
% 
%   F./G is the same as the command rdivide(F,G)
%
% See also LDIVIDE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% If either f or g are empty then return an empty SEPARABLEAPPROX object.
if ( isempty( f ) )
    h = f;
    return;
elseif ( isempty( g ) )
    h = g;
    return 
end

if ( isa(f, 'separableApprox') && isa(g, 'separableApprox') )    % SEPARABLEAPPROX ./ SEPARABLEAPPROX
    if ( ~domainCheck(f, g))
       error('CHEBFUN:SEPARABLEAPPROX:rdivide:domains', 'Domains inconsistent.') 
    end
    
    h = compose( f, @rdivide, g );
    
elseif ( isa(f, 'separableApprox') && isa(g, 'double') )  % SEPARABLEAPPROX ./ double 
    if ( g == 0 )
        error('CHEBFUN:SEPARABLEAPPROX:rdivide:divByZero', ...
            'Division by zero or near zero.')
    end
    h = f.* ( 1 / g ) ;
        
elseif ( isa(f, 'double') && isa(g, 'separableApprox') )   
       [bol, wzero] = singleSignTest( g );  
       if ( ( bol == 1 ) && ( wzero == 0 ) )
           h = compose( f, @rdivide, g );
       else
          error('CHEBFUN:SEPARABLEAPPROX:rdivide:zero', ...
              'Attempting to invert a SEPARABLEAPPROX with a root.'); 
       end
       
else
    error('CHEBFUN:CHEBFUN2:rdivide:badInputs', 'Unrecognised operation.');
    
end

end
