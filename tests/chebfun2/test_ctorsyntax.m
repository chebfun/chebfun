function pass = test_ctorsyntax( pref )
% This tests the Chebfun2 constructor for different syntax.
% Alex Townsend, March 2013.

pass(1) = 1;
f = @(x,y) cos(x) + sin(x.*y);  % simple function.
fstr = 'cos(x) + sin(x.*y)'; % string version.

try
    % % Adaptive calls % %
    % Operator way
    chebfun2(f);
    % String
    chebfun2(fstr);
    % With domain.
    chebfun2(f,[-1 1 1 2]);
    % Split domain syntax
    chebfun2(f,[-1,0],[2,3]);
    % Operator only in one variable.
    chebfun2(@(x,y) x);
    % Number pf 
catch
    pass(1) = 0 ;
end

% Test difference ways to make Chebfun2 nonadaptive. 
domain = [-1 1 -1 1]; 
f = @(x,y) sin(x.*y);
[xx, yy] = chebpts2( 10, 15, domain ); 
g = chebfun2( f(xx,yy), domain ); 
[m, n] = length(g);
pass(2) = ( m == 10 ) && (n == 15) ; 

f = @(x,y) sin(x.*y);
g = chebfun2( f, [10 15]);
[m,n] = length(g);
pass(3) = ( m == 10 ) && (n == 15) ; 

f = @(x,y) sin(x.*y);
g = chebfun2( f, [10 15], [-2 1 -2 1]);
[m,n] = length(g);
pass(4) = ( m == 10 ) && (n == 15) ;
pass(5) = any( g.domain == [-2 1 -2 1] ); 


end