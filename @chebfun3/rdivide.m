function h = rdivide(f,g)
%./   Pointwise CHEBFUN3 right divide.
%   F./G if F is a CHEBFUN3 and G is a double this returns (1/G)*F
%
%   F./G if F is a double and G is a v this returns F/G, but this does
%   not work if G becomes numerically close to zero.
%
%   F./G we do not allow F and G to both be CHEBFUN2 object.
% 
%   F./G is the same as the command rdivide(F,G)
%
% See also LDIVIDE.

% Empty check: 
if ( isempty( f ) || isempty( g ) ) 
    h = chebfun3();
    return
end

if ( isa(f, 'chebfun3') && isa(g, 'chebfun3') )    % CHEBFUN3 ./ CHEBFUN3
    if ( ~domainCheck(f, g))
       error('CHEBFUN:CHEBFUN3:rdivide:domains', 'Domains inconsistent.') 
    end
    h = chebfun3(@(x,y,z) feval(f,x,y,z)./feval(g,x,y,z), f.domain, ...
        'fiberDim', 3 ); 
    
elseif ( isa(f, 'chebfun3') && isa(g, 'double') )  % CHEBFUN3 ./ double 
    if ( g == 0 )
        error('CHEBFUN:CHEBFUN3:rdivide:divByZero', ...
            'Division by zero or near zero.')
    end
    h = f .* (1 / g) ;
        
elseif ( isa(f, 'double') && isa(g, 'chebfun3') )   
       [bool, wzero] = singleSignTest(g);
       if ( ( bool == 1 ) && ( wzero == 0 ) )
           h = chebfun3(@(x,y,z) f ./ feval(g, x, y, z), g.domain, ...
               'fiberDim', 3);
       else
          error('CHEBFUN:CHEBFUN3:rdivide:zero', ...
              'Attempting to invert a CHEBFUN3 with a root.'); 
       end
       
elseif ( isa(f,'chebfun3') && isa(g,'chebfun3v') )
    % TODO: RDIVIDE on the components: 
    
else
    error('CHEBFUN:CHEBFUN3:rdivide:badInputs', 'Unrecognised operation.');
    
end

end