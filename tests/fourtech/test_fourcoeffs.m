% Test file for fourtech/fourcoeffs.m

function pass = test_fourcoeffs(pref)

% Get preferences.
if ( nargin < 1 )
    pref = fourtech.techPref();
end

testclass = fourtech();

%%
% Check a few simple examples.

f = testclass.make(@(x) zeros(size(x)), [], pref);
p = fourcoeffs(f);
pass(1) = (norm(p, inf) <= 10*f.vscale.*f.epslevel);

f = testclass.make(@(x) 3*ones(size(x)), [], pref);
p = fourcoeffs(f);
pass(2) = (norm(p - 3, inf) < 10*f.vscale.*f.epslevel);

f = testclass.make(@(x) 1+cos(pi*x), [], pref);
p = fourcoeffs(f);
pass(3) = (norm(p - [0.5 1 0.5]', inf) < 10*f.vscale.*f.epslevel);
p = fourcoeffs(f,5);
pass(4) = (norm(p - [0 0.5 1 0.5 0]', inf) < 10*f.vscale.*f.epslevel);
p = fourcoeffs(f,1);
pass(5) = (norm(p - 1, inf) < 10*f.vscale.*f.epslevel);

f = testclass.make(@(x) 1 + exp(2*1i*pi*x) + exp(-1i*pi*x), [], pref);
p = fourcoeffs(f);
pass(6) = (norm(p - [0 1 1 0 1]', inf) ...
    < 10*f.vscale.*f.epslevel);
p = fourcoeffs(f,9);
pass(7) = (norm(p - [0 0 0 1 1 0 1 0 0]', inf) ...
    < 10*f.vscale.*f.epslevel);
p = fourcoeffs(f,3);
pass(8) = (norm(p - [1 1 0]', inf) ...
    < 10*f.vscale.*f.epslevel);

%%
% Verify operation for array-valued fourtech objects.

f = testclass.make(@(x) [3*ones(size(x)), 1+cos(pi*x), ... 
    1 + exp(2*1i*pi*x) + exp(-1i*pi*x)], [], pref);
p = fourcoeffs(f);
p_exact = [0 0   0;...
           0 0.5 1;...   
           3 1   1;...
           0 0.5 0;...
           0 0   1];
pass(9) = (norm(p(:) - p_exact(:), inf) < 10*max(f.vscale.*f.epslevel));

p = fourcoeffs(f,7);
p_exact = [0 0   0;...
           0 0   0;...
           0 0.5 1;...   
           3 1   1;...
           0 0.5 0;...
           0 0   1;...
           0 0   0];
pass(10) = (norm(p(:) - p_exact(:), inf) < 10*max(f.vscale.*f.epslevel));

p = fourcoeffs(f,3);
p_exact = [0 0.5 1;...   
           3 1   1;...
           0 0.5 0];
pass(11) = (norm(p(:) - p_exact(:), inf) < 10*max(f.vscale.*f.epslevel));

p = fourcoeffs(f,0);
pass(12) = isempty(p);

end
