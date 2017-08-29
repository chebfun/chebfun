function pass = test_iszero() 
% Test diskfun iszero() command.

f = diskfun( @(x,y) 0 + 0*x ); 
pass(1) = iszero( f ); 

f = diskfun( @(x,y) cos(x) );
pass(2) = ~iszero( f ); 

end 