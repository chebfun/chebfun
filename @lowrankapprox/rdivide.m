function h = rdivide(f,g)
%./   Pointwise right divide of LOWRANKAPPROX objects.
%   F./G if F is a LOWRANKAPPROX and G is a double this returns (1/G)*F
%
%   F./G if F is a double and G is a v this returns F/G, but this does
%   not work if G becomes numerically close to zero.
%
%   F./G we do not allow F and G to both be LOWRANKAPPROX object.
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

if ( isa(f, 'lowrankapprox') && isa(g, 'lowrankapprox') )    % LOWRANKAPPROX ./ LOWRANKAPPROX
    if ( ~domainCheck(f, g))
       error('CHEBFUN:LOWRANKAPPROX:rdivide:domains', 'Domains inconsistent.') 
    end
    %h = chebfun2( @(x,y) feval(f,x,y)./feval(g,x,y) , f.domain ); 
    
elseif ( isa(f, 'lowrankapprox') && isa(g, 'double') )  % LOWRANKAPPROX ./ double 
    if ( g == 0 )
        error('CHEBFUN:LOWRANKAPPROX:rdivide:divByZero', ...
            'Division by zero or near zero.')
    end
    h = f.* ( 1 / g ) ;
        
elseif ( isa(f, 'double') && isa(g, 'lowrankapprox') )   
       [bol, wzero] = singleSignTest( g );  
       if ( ( bol == 1 ) && ( wzero == 0 ) )
           %h = chebfun2( @(x,y) f ./ feval(g, x, y) , g.domain );
           
       else
          error('CHEBFUN:LOWRANKAPPROX:rdivide:zero', ...
              'Attempting to invert a LOWRANKAPPROX with a root.'); 
       end
       
elseif ( isa(f,'lowrankapprox') && isa(g,'lowrankapproxv') )
    % TODO: RDIVIDE on the components: 
    
else
    error('CHEBFUN:CHEBFUN2:rdivide:badInputs', 'Unrecognised operation.');
    
end

end
