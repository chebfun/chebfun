function pass = test_qr( pref ) 
% Test for QR decomposition of a chebfun2. 

if ( nargin == 0 ) 
    pref = chebfunpref; 
end 

tol = 100*pref.cheb2Prefs.chebfun2eps;
j = 1; 

% Decomposition on [-1,1,-1,1]: 
x = chebpts(100);
f = chebfun2( @(x,y) 1 ./ ( 10 + x.^2 + y.^2 ) );
[Q, R, E] = qr( f );
g = Q * R;

% Accurate decomposition: 
pass(j) = norm( fevalm(f, x, x) - fevalm(g, x, x) ) < sqrt( tol ); j = j + 1; 
pass(j) = norm( Q' * Q - eye( size(Q, 2) ) ) < tol; j = j + 1; 
pass(j) = norm( tril( R( :, E ) , -1) ) < tol; j = j + 1; 

%%

% Try the same thing on different domain: 
x = chebpts(100, [-2.1 4.3]);
y = chebpts(100, [-1 .3]);

f = chebfun2( @(x,y) 1 ./ ( 10 + x.^2 + y.^2 ) , [-2.1 4.3 -1 .3]);
[Q, R, E] = qr( f );
g = Q * R; 

% Accurate decomposition: 
pass(j) = norm( fevalm(f, x, y) - fevalm(g, x, y) ) < sqrt( tol ); j = j + 1; 
pass(j) = norm( Q' * Q - eye( size(Q, 2) ) ) < tol; j = j + 1; 
pass(j) = norm( tril( R( :, E ) , -1) ) < tol; j = j + 1;  


% Different syntax for the same thing: 
Q1 = qr( f ); 
Q2 = qr( f, 0 ); 
[Q3, R3] = qr( f ); 
[Q4, R4] = qr( f, 0); 
[Q5, R5, E1] = qr( f, 'vector' ); 
[Q6, R6, E2] = qr( f ); 

pass(j) = norm( Q1 - Q2 ) < tol; j = j + 1; 
pass(j) = norm( Q2 - Q3 ) < tol; j = j + 1; 
pass(j) = norm( Q3 - Q4 ) < tol; j = j + 1; 
pass(j) = norm( Q4 - Q5 ) < tol; j = j + 1; 
pass(j) = norm( Q5 - Q6 ) < tol; j = j + 1; 
pass(j) = norm( R3 - R4 ) < tol; j = j + 1; 
pass(j) = norm( R4 - R5 ) < tol; j = j + 1; 
pass(j) = norm( R5 - R6 ) < tol; j = j + 1; 
pass(j) = norm( E1 - E2 ) < tol; j = j + 1;

%% Make sure strictly real functions have a strictly real QR: 
f = chebfun2(@(x,y) sin(exp(x.*y)));
[Q,R] = qr(f);
pass(j) = norm( imag( Q ) ) == 0; j = j + 1; 
pass(j) = norm( imag( R ) ) == 0; j = j + 1; 
end