function pass = test_subsassgn(pref)

%% Setup.
if ( nargin == 0 )
    pref = chebfunpref();
end

%% See #2280

V = chebmatrix({1,2});

V(1) = pi;
pass(1) = V{1} == pi;

V{1} = -pi;
pass(2) = V{1} == -pi;

end
