function pass = test_inherited( ) 
% Test inherited commands, abs, cos, cosh, conj, etc. 

tol = 1e2*chebfunpref().cheb2Prefs.chebfun2eps; 

% abs
f = diskfun( @(x,y) -1+0*x);
h = diskfun( @(x,y) 1+0*x);
g = abs( f );
pass(1) = norm( g - h ) < tol;

% cos 
f = diskfun( @(x,y) cos(x));
h = diskfun( @(x,y) cos(cos(x)));
g = cos( f );
pass(2) = norm( g - h ) < tol;

% cosh
f = diskfun( @(x,y) cos(x));
h = diskfun( @(x,y) cosh(cos(x)));
g = cosh( f );
pass(3) = norm( g - h ) < tol;

% conj
f = diskfun( @(x,y) x);
g = conj( f );
pass(4) = norm( f - g ) < tol;


% exp
f = diskfun( @(x,y) cos(x)+sin(y));
h = diskfun( @(x,y) exp(cos(x)+sin(y)));
g = exp( f );
pass(5) = norm( g - h ) < tol;


% imag 
f = diskfun( @(x,y) x );
g = imag( f );
pass(6) = norm( g  ) < tol;

% isequal 
f = diskfun( @(x,y) cos(x));
g = diskfun( @(x,y) cos(-x) );
pass(7) = isequal( f, g ); 

f = diskfun( @(x,y) cos(x));
g = diskfun( @(x,y) cos(x) + 1 );
pass(8) = ~isequal( f, g ); 

% isreal 
f = diskfun( @(x,y) x);
pass(9) = isreal(f);

f = diskfun( @(x,y) cos(x.*y));
pass(10) = isreal(f); 

% iszero 
f = diskfun( @(x,y) 0*cos(x.*y));
pass(11) = iszero( f ); 
f = diskfun( zeros(10) ); 
pass(12) = iszero( f );

% length
f = diskfun( @(t,r) r.^2, 'polar');
[m, n] = length( f ); 
pass(13) = (m == 1);

% log 
f = diskfun( @(x,y) exp(x));
h = diskfun( @(x,y) x);
g = log( f );
pass(14) = norm( g - h ) < tol;


% sin 
f = diskfun( @(x,y) cos(x.*y));
h = diskfun( @(x,y) sin(cos(x.*y)));
g = sin( f );
pass(15) = norm( g - h ) < tol;

% sinh 
f = diskfun( @(x,y) cos(x.*y));
h = diskfun( @(x,y) sinh(cos(x.*y)));
g = sinh( f );
pass(16) = norm( g - h ) < tol;

% size 
pass(17) = (size( f,1 ) == inf && size( f,2 ) == inf);

% sqrt 
f = diskfun( @(x,y) cos(x.*y).^2);
h = diskfun( @(x,y) cos(x.*y));
g = sqrt( f );
pass(18) = norm( g - h ) < tol;

% tan 
f = diskfun( @(x,y) cos(x.*y));
h = diskfun( @(x,y) tan(cos(x.*y)));
g = tan( f );
pass(19) = norm( g - h ) < 10*tol;

% tand 
f = diskfun( @(x,y) cos(x.*y));
h = diskfun( @(x,y) tand(cos(x.*y)));
g = tand( f );
pass(20) = norm( g - h ) < tol;

% tanh
f = diskfun( @(x,y) cos(x.*y));
h = diskfun( @(x,y) tanh(cos(x.*y)));
g = tanh( f );
pass(21) = norm( g - h ) < tol;

% uminus 
f = diskfun( @(x,y) cos(x.*y));
h = diskfun( @(x,y) -(cos(x.*y)));
g = uminus( f );
pass(22) = norm( g - h ) < tol;

% uplus 
f = diskfun( @(x,y) cos(x.*y));
h = diskfun( @(x,y) +(cos(x.*y)));
g = uplus( f );
pass(23) = norm( g - h ) < tol;

end