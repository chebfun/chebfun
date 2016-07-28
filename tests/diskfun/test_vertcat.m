function pass = test_vertcat( )
% test spherefun vertcat() command

f = spherefun( @(x,y,z) cos(x) );
F = [ f ; f ; f ];   % check that this works
pass(1) = iszero(F(1) - f);

% Note that this does not work:
f = spherefun( @(x,y,z) cos(x) );
try
    F = [ f ; f ];
    pass(2) = 0;
catch
    pass(2) = 1;
end

end