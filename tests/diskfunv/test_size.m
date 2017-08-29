function pass = test_size( pref ) 
% Test size

% Check three components, which is all that is allowed right now
F = diskfunv(@(x,y) cos(x), @(x,y) sin(y)); 
pass(1) = all( size(F) == [2 inf inf]); 
pass(2) = ( size(F, 1) == 2 ); 
pass(3) = ( size(F, 2) == inf); 
pass(4) = ( size(F, 3) == inf);

%check transpose
F = F';
pass(5) = all( size(F) == [inf inf 2]); 
pass(6) = ( size(F, 1) == inf ); 
pass(7) = ( size(F, 2) == inf); 
pass(8) = ( size(F, 3) == 2);

end