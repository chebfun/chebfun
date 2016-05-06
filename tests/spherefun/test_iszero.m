function pass = test_iszero() 
% Test spherefun iszero() command.

f = spherefun( @(x,y,z) 0 + 0*x ); 
pass(1) = iszero( f ); 

f = spherefun( @(x,y,z) cos(x) );
pass(2) = ~iszero( f ); 

end 