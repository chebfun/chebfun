function pass = test_iszero( pref ) 

% Grab some preferences
if ( nargin == 0 )
    pref = chebfunpref();
end

% Test with vector 1-1
f = ballfun(ones(19,20,22));
F = ballfunv(f,f,f);
G = F-F;
pass(1) = iszero(G);
end
