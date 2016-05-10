function pass = test_roots( )
% Test the chebfun2/roots command.

tol = 1e-14;

% Here, we check that the roots command is working for separable functions. 

f = chebfun2( @(x,y) 1/2-y);
r = roots( f );
exact = chebfun( @(x) x + 1i/2 );
pass(1) = ( norm( r - exact ) < tol );

f = chebfun2( @(x,y) 1/2-x);
r = roots( f );
exact = chebfun( @(x) 1i*x + 1/2 );
pass(2) = ( norm( r - exact ) < tol );


f = chebfun2( @(x,y) x.*y);
r = roots( f );
exact = chebfun( @(x) [x, 1i*x]  );
pass(3) = ( norm( r - exact ) < tol );


f = chebfun2( @(x,y) cos(5*pi*x));
s = roots( squeeze( f ) );
r = roots( f );
exact = chebfun;
for j = 1:numel(s)
    exact = [ exact chebfun(@(x) s(j) + 1i*x )];
end
pass(4) = ( norm( r - exact ) < tol );

f = chebfun2( @(x,y) x.*(y-1/2), [-2 2 -3 5]);
r = roots( f );
exact = [chebfun( @(x) 2*x + 1i/2 ) chebfun(@(x) 1i*(4*(x+1)-3) )];
pass(5) = ( norm( r - exact ) < tol );

end