function pass = test_isequal( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end

% Test with function 1
f = ballfun(ones(21,20,22));
g = f+f-f;

pass(1) = (isequal(f,g) && isequal(g,f));
end
