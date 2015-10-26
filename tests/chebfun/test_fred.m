function pass = test_fred(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

K = @(u, v) exp(-(u-v).^2);     % A simple kernel.
f = chebfun(@sin, pref);        % A simple chebfun.
F = fred(K, f);                 % Call fred().

% Test against some V4 results:
pass(1) = abs(F(.5) - 0.293968048825243) < 1e1*vscale(F)*eps;
pass(2) = abs(norm(F) - 0.392002900508830) < 1e1*vscale(F)*eps;


% Test 3rd input argument. (Simply make sure we don't crash!)
try 
    F2 = fred(K, f, 1);
    pass(3) = 1;
catch
    pass(3) = 0;
end

end
