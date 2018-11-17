function pass = test_iszero( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end

% Test with function 1-1
f = ballfun(ones(21,20,22));
g = f-f;

pass(1) = iszero(g);
end
