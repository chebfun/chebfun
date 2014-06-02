function pass = test_subsref( pref ) 
% Test CHEBOP2 SUBSREF 

N = chebop2(@(u) diff(u,2,1) + diff(u,2,2)); 

pass(1) = ( size( N(10,10), 1) - 100 ) == 0; 
pass(2) = ( size( N(10,10), 2) - 100 ) == 0;

m = 5; n = 10; 
pass(3) = ( size( N(m,n), 1) - m*n ) == 0;
pass(4) = ( size( N(m,n), 2) - m*n ) == 0;

m = 5; n = 10; 
pass(5) = ( size( N(m), 1) - m^2 ) == 0;
pass(6) = ( norm( N(n,n) - N(n) ) == 0 ); 


% GET a property: 
pass(7) = ( norm( N.coeffs - [0 0 1; 0 0 0; 1 0 0]) == 0 ); 


end
