function pass = test_volt(pref)

if ( nargin == 0 )
    pref = chebpref();
end

% A simple kernel:
K = @(u, v) exp(-(u-v).^2);
% Asimple CHEBFUN:
f = chebfun(@sin, pref);
% Call VOLT
F = volt(K, f);

% Test against some V4 results:
pass(1) = abs(F(.5) - (-0.013808570536509)) < vscale(F)*epslevel(F);
pass(2) = abs(norm(F) - 0.334612395278957) < vscale(F)*epslevel(F);

% Test 3rd input argument. (Simply make sure we don't crash!)
try 
    F2 = volt(K, f, 1);
    pass(3) = 1;
catch
    pass(3) = 0;
end

end
