function pass = test_size( pref ) 
% Test size

% Check three components, which is all that is allowed right now
F = spherefunv(@(x,y,z) cos(x), @(x,y,z) sin(y), @(x,y,z) cos(1+z)); 
pass(1) = all( size(F) == [3 inf inf]); 
pass(2) = ( size(F, 1) == 3 ); 
pass(3) = ( size(F, 2) == inf); 
pass(4) = ( size(F, 3) == inf); 

end