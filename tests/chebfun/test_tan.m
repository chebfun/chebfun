function pass = test_tan(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% Bug from #568.
f = chebfun('tan(x)', pi*((-5/2):(5/2)), pref, 'exps', -ones(1, 6));
pass = length(f) < 2^16;

end