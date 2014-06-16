function pass = test_volt(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

K = @(u, v) exp(-(u-v).^2);     % A simple kernel.
f = chebfun(@sin, pref);        % A simple chebfun.
F = volt(K, f);                 % Call volt().

% Test against some V4 results:
pass(1) = abs(F(.5) - (-0.013808570536509)) < 10*vscale(F)*epslevel(F);
pass(2) = abs(norm(F) - 0.334612395278957) < vscale(F)*epslevel(F);

% Test 3rd input argument. (Simply make sure we don't crash!)
try 
    F2 = volt(K, f, 1);
    pass(3) = 1;
catch
    pass(3) = 0;
end

end
