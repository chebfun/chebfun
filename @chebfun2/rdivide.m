function h = rdivide(f,g)
%./   Pointwise chebfun2 right divide.
%
% F./G if F is a chebfun2 and G is a double this returns (1/G)*F
%
% F./G if F is a double and G is a chebfun2 this returns F/G, but this
% does not work if G becomes numerically close to zero.
%
% F./G we do not allow F and G to both be chebfun2 object.
% 
% F./G is the same as the command rdivide(F,G)
%
% See also LDIVIDE.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check: 
if ( ( isempty( f ) ) || ( isempty( g ) ) ) 
    h = chebfun2;
    return;
end

if ( isa(f,'chebfun2') && isa(g,'chebfun2') )    % chebfun2 ./ chebfun2 
    if ( ~chebfun2.domainCheck(f, g))
       error('CHEBFUN2:RDIVIDE:DOMAINS','Domains inconsistent.') 
    end
    h = chebfun2( @(x,y) feval(f,x,y)./feval(g,x,y) , f.domain ); 
elseif ( isa(f,'chebfun2') && isa(g,'double') )  % chebfun2 ./ double 
        if ( g == 0 )
            error('CHEBFUN2:rdivide:DivisionByZero','Division by zero or near zero.')
        end
        h = f.* ( 1 / g ) ;
elseif ( isa(f,'double') && isa(g,'chebfun2') )   
       [bol, wzero] = chebfun2.singleSignTest( g );  
       if ( ( bol == 1 ) && ( wzero == 0 ) )
           h = chebfun2( @(x,y) f ./ feval(g, x, y) , g.domain );
       else
          error('CHEBFUN2:RDIVIDE:ZERO','Attempting to invert a chebfun2 with a root.'); 
       end
elseif ( isa(f,'chebfun2') && isa(g,'chebfun2v') )
    % Do RDIVIDE on the components: 
    
    % TODO: 
else
    error('CHEBFUN2:rdivide:Inputs','Unrecognised operation.');
end


end