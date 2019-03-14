function pass = test_iszero( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end

% Test with function 1-1
f = ballfun(ones(21,20,22));
g = f-f;
pass(1) = iszero(g);

% Test with function 10^-20
f = ballfun(10^(-20));
pass(2) = (iszero(f) == 0);

% Test with function 0
f = ballfun(0);
pass(3) = iszero(f);
end
