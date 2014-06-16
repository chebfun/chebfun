function pass = test_size( pref ) 
% Test SIZE

if ( nargin == 0 ) 
    pref = chebfunpref; 
end

% Check two components: 
F = chebfun2v(@(x,y) cos(x), @(x,y) sin(y));
pass(1) = all( size(F) == [2 inf inf]); 
pass(2) = ( size(F, 1) == 2 ); 
pass(3) = ( size(F, 2) == inf); 
pass(4) = ( size(F, 3) == inf); 

% Check three components: 
F = chebfun2v(@(x,y) cos(x), @(x,y) sin(y), @(x,y) cos(x)); 
pass(5) = all( size(F) == [3 inf inf]); 
pass(6) = ( size(F, 1) == 3 ); 
pass(7) = ( size(F, 2) == inf); 
pass(8) = ( size(F, 3) == inf); 

end