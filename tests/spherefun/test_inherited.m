function pass = test_inherited( ) 
% Test inherited commands, abs, cos, cosh, conj, etc. 

tol = 1e2*chebfunpref().cheb2Prefs.chebfun2eps; 

% abs
f = spherefun( @(x,y,z) -1+0*x);
h = spherefun( @(x,y,z) 1+0*x);
g = abs( f );
pass(1) = norm( g - h ) < tol;

% cos 
f = spherefun( @(x,y,z) cos(x));
h = spherefun( @(x,y,z) cos(cos(x)));
g = cos( f );
pass(2) = norm( g - h ) < tol;

% cosh
f = spherefun( @(x,y,z) cos(x));
h = spherefun( @(x,y,z) cosh(cos(x)));
g = cosh( f );
pass(3) = norm( g - h ) < tol;

% conj
f = spherefun( @(x,y,z) x);
g = conj( f );
pass(4) = norm( f - g ) < tol;

% ctranspose
f = spherefun( @(x,y,z) y);
g = ctranspose( f );
pass(5) = norm( f - g ) < tol;

% diag 
f = spherefun( @(lam,theta) sin(theta).*cos(lam));
h = chebfun( @(x) sin(x).*cos(x), [0 pi]);
g = diag( f );
pass(6) = norm( g - h ) < tol;

% exp
f = spherefun( @(x,y,z) cos(x)+sin(y));
h = spherefun( @(x,y,z) exp(cos(x)+sin(y)));
g = exp( f );
pass(7) = norm( g - h ) < tol;

% flipdim, fliplr, flipud
f = spherefun( @(lam, theta) sin(theta).*cos(lam));
g1 = spherefun( @(lam, theta) sin(theta).*cos(-lam));
g2 = spherefun( @(lam, theta) sin(-theta).*cos(lam));
h1 = fliplr( f ); 
h2 = flipud( f ); 
h3 = flipdim( f, 1); 
h4 = flipdim( f, 2); 
pass(8) = ( norm( g1 - h1 ) < tol );
pass(9) = ( norm( g2 - h2 ) < tol );
pass(10) = ( norm( g2 - h3 ) + norm( g1 - h4 ) < tol ); 

% imag 
f = spherefun( @(x,y,z) x );
g = imag( f );
pass(11) = norm( g  ) < tol;

% isequal 
f = spherefun( @(x,y,z) cos(x));
g = spherefun( @(x,y,z) cos(-x) );
pass(12) = isequal( f, g ); 

f = spherefun( @(x,y,z) cos(x));
g = spherefun( @(x,y,z) cos(x) + 1 );
pass(13) = ~isequal( f, g ); 

% isreal 
f = spherefun( @(x,y,z) x);
pass(14) = isreal(f);

f = spherefun( @(x,y,z) cos(x.*y.*z));
pass(15) = isreal(f); 

% iszero 
f = spherefun( @(x,y,z) 0*cos(x.*y.*z));
pass(16) = iszero( f ); 
f = spherefun( zeros(10) ); 
pass(17) = iszero( f );

% length
f = spherefun( @(x,y,z) cos(z));
[m, n] = length( f ); 
pass(18) = (m == 1);

% log 
f = spherefun( @(x,y,z) exp(z));
h = spherefun( @(x,y,z) z);
g = log( f );
pass(19) = norm( g - h ) < tol;

% rank

% sin 
f = spherefun( @(x,y,z) cos(x.*y));
h = spherefun( @(x,y,z) sin(cos(x.*y)));
g = sin( f );
pass(20) = norm( g - h ) < tol;

% sinh 
f = spherefun( @(x,y,z) cos(x.*y));
h = spherefun( @(x,y,z) sinh(cos(x.*y)));
g = sinh( f );
pass(21) = norm( g - h ) < tol;

% size 
pass(22) = (size( f,1 ) == inf && size( f,2 ) == inf);

% sqrt 
f = spherefun( @(x,y,z) cos(x.*y).^2);
h = spherefun( @(x,y,z) cos(x.*y));
g = sqrt( f );
pass(23) = norm( g - h ) < tol;

% tan 
f = spherefun( @(x,y,z) cos(x.*y));
h = spherefun( @(x,y,z) tan(cos(x.*y)));
g = tan( f );
pass(24) = norm( g - h ) < 10*tol;

% tand 
f = spherefun( @(x,y,z) cos(x.*y));
h = spherefun( @(x,y,z) tand(cos(x.*y)));
g = tand( f );
pass(25) = norm( g - h ) < tol;

% tanh
f = spherefun( @(x,y,z) cos(x.*y));
h = spherefun( @(x,y,z) tanh(cos(x.*y)));
g = tanh( f );
pass(26) = norm( g - h ) < tol;

% uminus 
f = spherefun( @(x,y,z) cos(x.*y));
h = spherefun( @(x,y,z) -(cos(x.*y)));
g = uminus( f );
pass(27) = norm( g - h ) < tol;

% uplus 
f = spherefun( @(x,y,z) cos(x.*y));
h = spherefun( @(x,y,z) +(cos(x.*y)));
g = uplus( f );
pass(28) = norm( g - h ) < tol;

end