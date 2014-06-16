function pass = test_ctorsyntax( pref )
% This tests the Chebfun2 constructor for different syntax.
% Alex Townsend, March 2013.

pass = 1;
f = @(x,y) cos(x) + sin(x.*y);  % simple function.
fstr = 'cos(x) + sin(x.*y)'; % string version.

try
    % % Adaptive calls % %
    % Operator way
    chebfun2(f);
    % String
    chebfun2(fstr);
    % With domain.
    chebfun2(f,[-1 1 1 2]);
    % Split domain syntax
    chebfun2(f,[-1,0],[2,3]);
    % Operator only in one variable.
    chebfun2(@(x,y) x);
catch
    pass = 0 ;
end

end