function h = rdivide(f,g)
%./   Pointwise CHEBFUN2 right divide.
%   F./G if F is a CHEBFUN2 and G is a double this returns (1/G)*F
%
%   F./G if F is a double and G is a v this returns F/G, but this does
%   not work if G becomes numerically close to zero.
%
%   F./G we do not allow F and G to both be CHEBFUN2 object.
% 
%   F./G is the same as the command rdivide(F,G)
%
% See also LDIVIDE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty( f ) || isempty( g ) ) 
    h = chebfun2();
    return
end

if ( isa(f, 'chebfun2') && isa(g, 'chebfun2') )    % CHEBFUN2 ./ CHEBFUN2 
    if ( ~domainCheck(f, g))
       error('CHEBFUN:CHEBFUN2:rdivide:domains', 'Domains inconsistent.') 
    end
    h = chebfun2( @(x,y) feval(f,x,y)./feval(g,x,y) , f.domain ); 
    
elseif ( isa(f, 'chebfun2') && isa(g, 'double') )  % CHEBFUN2 ./ double 
    if ( g == 0 )
        error('CHEBFUN:CHEBFUN2:rdivide:divByZero', ...
            'Division by zero or near zero.')
    end
    h = f.* ( 1 / g ) ;
        
elseif ( isa(f, 'double') && isa(g, 'chebfun2') )   
       [bol, wzero] = singleSignTest( g );  
       if ( ( bol == 1 ) && ( wzero == 0 ) )
           h = chebfun2( @(x,y) f ./ feval(g, x, y) , g.domain );
       else
          error('CHEBFUN:CHEBFUN2:rdivide:zero', ...
              'Attempting to invert a CHEBFUN2 with a root.'); 
       end
       
elseif ( isa(f,'chebfun2') && isa(g,'chebfun2v') )
    % TODO: RDIVIDE on the components: 
    
else
    error('CHEBFUN:CHEBFUN2:rdivide:badInputs', 'Unrecognised operation.');
    
end

end
