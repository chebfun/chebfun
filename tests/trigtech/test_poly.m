% Test file for trigtech/poly.m

function pass = test_poly(pref)

% Get preferences.
if ( nargin < 1 )
    pref = trigtech.techPref();
end

testclass = trigtech();

%%
% Check a few simple examples.

f = testclass.make(@(x) zeros(size(x)), [], pref);
p = poly(f);
pass(1) = (norm(p, inf) <= 10*vscale(f)*eps);

f = testclass.make(@(x) 3*ones(size(x)), [], pref);
p = poly(f);
pass(2) = (norm(p - 3, inf) < 10*vscale(f)*eps);

f = testclass.make(@(x) 1+cos(pi*x), [], pref);
p = poly(f);
pass(3) = (norm(p - [0.5 1 0.5], inf) < 10*vscale(f)*eps);

f = testclass.make(@(x) 1 + exp(2*1i*pi*x) + exp(-1i*pi*x), [], pref);
p = poly(f);
pass(4) = (norm(p - [0 1 1 0 1], inf) ...
    < 10*vscale(f)*eps);

%%
% Verify operation for array-valued chebtech objects.

f = testclass.make(@(x) [3*ones(size(x)), 1+cos(pi*x), ... 
    1 + exp(2*1i*pi*x) + exp(-1i*pi*x)], [], pref);
p = poly(f);
p_exact = [0 0   0;...
           0 0.5 1;...   
           3 1   1;...
           0 0.5 0;...
           0 0   1].';
pass(5) = (norm(p(:) - p_exact(:), inf) < 10*max(vscale(f)*eps));

end
