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

if isempty(f) || isempty(g)  % isempty check. 
    h = chebfun2;
    return;
end


if ( isa(f,'chebfun2') && isa(g,'chebfun2') )
%% chebfun2 ./ chebfun2 
    rectf = f.corners; rectg = g.corners; 
    if ( norm(rectf - rectg,inf) >0 ) 
        error('CHEBFUN2:rdivide:domain','Inconsistent domains')
    end
    h = chebfun2( @(x,y) feval(f,x,y)./feval(g,x,y) , rectf ); 
elseif ( isa(f,'chebfun2') && isa(g,'double') )
%% chebfun2 ./ double 
        if ( abs(g) < eps )
            error('CHEBFUN2:rdivide:DivisionByZero','Division by zero or near zero.')
        end
        h = f.*(1/g);
elseif ( isa(f,'double') && isa(g,'chebfun2') )
       [bol wzero] = singlesigntest(g);   % only attempt if g is of single sign. 
       if bol == 1 && wzero == 0
           h = chebfun2( @(x,y) f./feval(g,x,y) ,g.corners );
       else
          error('CHEBFUN2:RDIVIDE:ZERO','Attempting to invert a chebfun2 with a zero.'); 
       end
elseif ( isa(f,'chebfun2') && isa(g,'chebfun2v') )
    h=g; 
    h.xcheb = rdivide(f,g.xcheb);
    h.ycheb = rdivide(f,g.ycheb);
    h.zcheb = rdivide(f,g.zcheb);
else
    error('CHEBFUN2:rdivide:Inputs','Unrecognised operation.');
end


end