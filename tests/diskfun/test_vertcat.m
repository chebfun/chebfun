function pass = test_vertcat( )
% test diskfun vertcat() command

f = diskfun( @(x,y) cos(x) );
F = [ f ; f ];   % check that this works
pass(1) = iszero(F(1) - f);

% Note that this does not work:
f = diskfun( @(x,y) cos(x) );
try
    F = [ f ; f; f ];
    pass(2) = 0;
catch
    pass(2) = 1;
end

end