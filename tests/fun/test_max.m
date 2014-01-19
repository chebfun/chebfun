% Test file for fun/max.m

function pass = test_max(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebpref();
end

% Set a domain for BNDFUN.
dom = [-2 7];
    
%% 
% Spot-check the extrema for a few BNDFUN.
pass(1) = test_spotcheck_max(@(x) sin(10*x), dom, 1, pref);
pass(2) = test_spotcheck_max(@airy, dom, 0.535656656015700, pref);
pass(3) = test_spotcheck_max(@(x) -1./(1 + x.^2), dom, -.02, pref);
pass(4) = test_spotcheck_max(@(x) (x/10).^3.*cosh(x/10), ...
    dom, 0.7^3*cosh(0.7), pref);

%%
% Check operation for array-valued BNDFUN inputs.
fun_op = @(x) [sin(10*x) airy(x) (x/10).^3.*cosh(x/10)];
f = bndfun(fun_op, dom, [], [], pref);
[y, x] = max(f);
exact_max = [1 0.535656656015700 0.7^3*cosh(0.7)];
fx = [sin(10*x(1)) airy(x(2)) (x(3)/10).^3.*cosh(x(3)/10)];
tol = 10*get(f, 'vscale').*get(f, 'epslevel');
pass(5) = (all(abs(y - exact_max) < tol) && ...
    all(abs(fx - exact_max) < tol));

%%
% Test for complex-valued BNDFUN.
pass(6) = test_spotcheck_max(@(x) (x/2).*(exp(1i*(x/2))+1i*sin(x/2)), ...
    dom, -3.277598405517787 - 2.455482593827339i, pref);

fun_op = @(x) [((x-2).^2/4+1).*exp(1i*(x/2)) ... 
    -((x+1).^2/4+1).*exp(1i*(x/2))];
f = bndfun(fun_op, dom, [], [], pref);
[y, x] = max(f);
exact_max = [-6.789310982858273-2.543178400749744i 15.919763683943538+5.963314870723537i];
fx = [((x(1)-2).^2/4+1).*exp(1i*(x(1)/2)) ... 
    -((x(2)+1).^2/4+1).*exp(1i*(x(2)/2))];
tol = get(f, 'vscale').*get(f, 'epslevel');
pass(7) = (all(abs(y - exact_max) < tol) && ...
    all(abs(fx - exact_max) < tol));
      
%% 
% [TODO]: Run a few tests for UNBNDFUN.
end

%%
% Spot-check the results for a given BNDFUN.
function result = test_spotcheck_max(fun_op, dom, exact_max, pref)

f = bndfun(fun_op, dom, [], [], pref);
[y, x] = max(f);
fx = fun_op(x);
tol = 10*get(f, 'vscale').*get(f, 'epslevel');
result = (all(abs(y - exact_max) < tol) && ...
          all(abs(fx - exact_max) < tol));

end