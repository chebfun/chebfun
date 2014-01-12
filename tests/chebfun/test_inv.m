function pass = test_inv(pref)
% This test constructs two CHEBFUN objects and uses INV() to invert them. It
% checks the inverse calculated is accurate.

% Taken from Chebfun v4 test, invtest.m, by Nick Hale  07/06/2009.

if ( nargin == 0 ) 
    pref = chebpref();
end

funcList = {@inv, @inv2};

for k = 1:2
    myinv = funcList{k};
    
    x = chebfun('x');
    f = sin(x);
    g = chebfun(@(x) asin(x), [sin(-1), sin(1)]);
    f_inv = myinv(f, pref);
    tol = 100*epslevel(f_inv).*vscale(f_inv);
    pass(k,1) = norm(g - f_inv,inf) < tol;

    pass(k,2) = norm(f - myinv(f_inv), inf) < tol;

    x = chebfun('x');
    f = chebfun(@(x) sausagemap(x));
    f_inv = myinv(f, pref);
    tol = 100*epslevel(f_inv).*vscale(f_inv);
    pass(k,3) = norm(f(f_inv) - x, inf) + norm(f_inv(f) - x, inf) < tol;
end
    

end

function g = sausagemap(s,d)
if ( nargin < 2 )
    d = 9; % This can be adjusted
end 
c = zeros(1,d+1);
c(d:-2:1) = [1 cumprod(1:2:d-2)./cumprod(2:2:d-1)]./(1:2:d);
c = c/sum(c); g = polyval(c,s);
cp = c(1:d).*(d:-1:1);
end
