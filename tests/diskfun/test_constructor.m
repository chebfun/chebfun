function pass = test_constructor( ) 
% Test the diskfun constructor 

% Get tolerance: 
tol = 2e3*chebfunpref().cheb2Prefs.chebfun2eps;

f = @(t,r) (r.*cos(t)).^2 + (r.*sin(t)).^2 ; 
g = diskfun( f, 'polar' );
pass(1) = ( SampleError( f,g ) < tol ); 

f = @(x,y) exp(-cos(pi*(x+y))); 
g = diskfun( f );
f = redefine_function_handle( f);
pass(2) = ( SampleError( f, g ) < tol ); 

f = @(t,r) cos(pi*r.*cos(t))+sin(5*(r.*sin(t)))-1; 
g = diskfun( f, 'polar' );
pass(3) = ( SampleError( f, g ) < tol ); 

f = @(x,y) exp(y); % is it ok to throw in 'cart', even though not needed?
g = diskfun( f , 'cart');
f = redefine_function_handle( f);
pass(4) = ( SampleError( f, g ) < tol ); 

f = @(x,y) 1-exp(x); 
g = diskfun( f );
f = redefine_function_handle( f);
pass(5) = ( SampleError( f, g ) < tol ); 

f = @(x,y) exp(-x)+exp(-y) ; 
g = diskfun( f );
f = redefine_function_handle( f );
pass(6) = ( SampleError( f, g ) < tol ); 

f = @(x,y) cos(2*x.*y); 
g = diskfun( f );
f = redefine_function_handle( f);
pass(7) = ( SampleError( f, g ) < tol ); 

f = @(x,y) sin(x.*y);
g = diskfun( f );
f = redefine_function_handle( f );
pass(8) = ( SampleError( f, g ) < tol ); 

f = @(x,y) sin(11*pi*x)-sin(3*pi*x) + sin(y);
g = diskfun( f );
f = redefine_function_handle( f );
pass(9) = ( SampleError( f, g ) < tol ); 

f = @(x,y) sech(cos(9*y)+sin(7*(x-1))); %rank 102
g = diskfun( f);
f = redefine_function_handle( f);
pass(10) = ( SampleError( f, g ) < 2e1*tol );


f = @(x,y) 0*x;
g = diskfun(f);
h = diskfun(f, 'polar');
pass(11) = ( norm(g, inf) == 0 ); 
pass(12) = (norm(h, inf) ==0);

% Test the vectorize flag is working: 
f = diskfun(@(x,y) cos(x));
g = diskfun(@(x,y) cos(x), 'vectorize');
pass(13) = ( norm(f - g) < tol );

f = diskfun(@(t,r) r.*sin(t), 'polar');
g = diskfun(@(t,r) r*sin(t), 'polar', 'vectorize');
pass(14) = ( norm(f - g) < tol );

f = diskfun(@(x,y) x.*y);
g = diskfun(@(x,y) x*y, 'vectorize');
pass(15) = ( norm(f-g) < tol ); 

f = diskfun(@(x,y) 1);
g = diskfun(@(x,y) 1, 'vectorize');
pass(16) = ( norm(f - g) < tol );

% Test construction from samples
f = diskfun(@(x,y) 1 + x.*sin((x-.5).*y));
[m,n] = length(f);
F = sample(f,m+mod(m,2),n);
g = diskfun(F);
pass(17) = ( norm(f - g) < tol );

f = diskfun(@(x,y) 1 + 0*x);
F = ones(1, 1);
g = diskfun(F);
pass(18) = ( norm(f - g) < tol );

% Test construction from coefficients.
f = diskfun(@(x,y) exp(-10*((x-.1).^2  + y.^2)));
C = coeffs2(f);
g = diskfun(C,'coeffs');
pass(19) = ( norm(f - g) < tol );

% Test fixed rank construction
ff = @(x,y) exp(-10*((x-.1).^2  + y.^2));
f = diskfun(ff,5);
pass(20) = ( rank(f) == 5 );
f = diskfun(ff,6);
pass(21) = ( rank(f) == 6 );
f = diskfun(ff);
g = diskfun(f,7);
pass(22) = ( rank(g) == 7 );
ff = @(t, r) exp(-10*((r.*cos(t)-.1).^2+(r.*sin(t)).^2));  %polar
f = diskfun(ff, 5, 'polar'); 
pass(23) = ( rank(f)== 5 ); 

% Test zero rank construction gives zero.
g = diskfun(f,0);
pass(24) = ( rank(g) == 0 );
pass(25) = norm(g) < tol;

try
    f = diskfun(ff,-1);
    pass(26) = false;
catch ME
    pass(26) = strcmp(ME.identifier,'CHEBFUN:DISKFUN:constructor:parseInputs:domain3');
end

% Check the 'eps' flag works.
ff = @(x,y) exp(-10*((x-.1).^2  + y.^2));
f = diskfun(ff);
g = diskfun(ff,'eps', 1e-5);
pass(27) = rank(g) < rank(f);
[mf,nf] = length(f);
[mg,ng] = length(g);
pass(28) = ( (mg < mf) && (ng < nf) );

% Construction from a single value
f = diskfun(1);
pass(29) = norm(f-1,inf) == 0;

% Construction from a string with (x,y) variables
f = diskfun(@(x,y) exp(-10*((x-.1).^2  + y.^2)));
g = diskfun('exp(-10*((x-.1).^2  + y.^2))');
pass(30) = norm(f-g,inf) == 0;

% Construction from a string with (t,r) variables
f = diskfun(@(t,r) r.*cos(t), 'polar' );
g = diskfun('y.*cos(x)', 'polar'); %alphabetic order; can't use r and t
pass(31) = norm(f-g,inf) == 0;

try
    f = diskfun('x.*y.*z');
    pass(32) = false;
catch ME
    pass(32) = strcmp(ME.identifier,'CHEBFUN:DISKFUN:constructor:str2op:depvars');
end

% Construction from zeros matrix should maintain a zeros coefficient
% matrix.
f = diskfun(zeros(5,4));
[n,m] = length(f);
pass(33) = (m == 5) && (n == 4);

% Fixed length test.
m = 24; n = 10;
f = diskfun(@(x,y) exp(-10*((x-0.5/sqrt(2)).^2 + (y-0.5/sqrt(2)).^2)),[m n]);
[nf,mf] = length(f);
pass(34) = (mf == m+1) && (nf == 2*n);

end

function sample_error = SampleError(h, g)
m = 7; 
n = m-1;
[x, y] = getPoints(m, n);
[L2, T2] = meshgrid(x, y);
F = h(L2, T2);
approx = fevalm(g, x, y);
sample_error = norm(F(:) - approx(:), inf);
end

function [x, y] = getPoints(m, n)

x = trigpts(2*n, [-pi pi]);
y = chebpts(m);
y = y(ceil(m/2):end); 

end

function f = redefine_function_handle(f)
    % Wrap Cartesian f so it can be evaluated in polar coordinates
    
    f = @(th, r) diskfun.pol2cartf(f,th, r);
  

end