function pass = test_isempty() 
% Test diskfun isempty() command.

f = diskfun( ); 
pass(1) = isempty( f ); 

f = diskfun( @(x,y) cos(x) );
pass(2) = ~isempty( f ); 

end 