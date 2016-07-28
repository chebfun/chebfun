function pass = test_isempty() 
% Test spherefun isempty() command.

f = spherefun( ); 
pass(1) = isempty( f ); 

f = spherefun( @(x,y,z) cos(x) );
pass(2) = ~isempty( f ); 

end 