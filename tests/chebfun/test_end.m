function pass = test_end(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
end

f = chebfun(@(x) sin(pi*x));
out = f(end);
pass = false;

% [TODO]: This requires SUBSREF() to work.

end