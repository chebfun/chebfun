% Test file for POLY.

function pass = test_poly(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

%% Test POLY on [-1 1]:

x = 0;
f = poly(x, domain([-1 1]));
err = feval(f, x);
pass(1) = err < 100*epslevel(f)*vscale(f);

x = 0.1;
f = poly(x, domain([-1 1]));
err = feval(f, x);
pass(2) = length(f) == 2 && err < 100*epslevel(f)*vscale(f);

x = [-.1 0.1];
f = poly(x, domain([-1 1]));
err = norm(feval(f, x));
pass(3) = length(f) == 3 && err < 100*epslevel(f)*vscale(f);

x = [-.1 0.1 ; -.5 .9];
f = poly(x, domain([-1 1]));
err = norm([f(x(:,1),1) f(x(:,2),2)]);
pass(4) = length(f) == 3 && size(f, 2) == 2 && err < 100*epslevel(f)*vscale(f);

x = [-.1 0.1 ; -.5 .9];
f = poly(x, domain([-1 2]));
err = norm([f(x(:,1),1) f(x(:,2),2)]);
pass(5) = length(f) == 3 && size(f, 2) == 2 && err < 100*epslevel(f)*vscale(f);

end