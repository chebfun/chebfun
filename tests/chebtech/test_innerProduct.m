% Test file for chebtech/innerProduct.m

function pass = test_innerProduct(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebtech.techPref();
end

% Fixed arbitrary numbers to use as multiplicative constants.
alpha = -0.194758928283640 + 0.075474485412665i;
beta = -0.526634844879922 - 0.685484380523668i;

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end

    %%
    % Spot-check a few known results.
    
    f = testclass.make(@(x) sin(2*pi*x), [], pref);
    g = testclass.make(@(x) cos(2*pi*x), [], pref);
    tol_f = 10*eps*vscale(f);
    tol_g = 10*eps*vscale(g);
    pass(n, 1) = abs(innerProduct(f, g)) < max(tol_f, tol_g);
    
    g = testclass.make(@(x) cos(4*pi*x), [], pref);
    tol_g = 10*eps*vscale(g);
    pass(n, 2) = abs(innerProduct(f, g)) < max(tol_f, tol_g);
    
    f = testclass.make(@(x) exp(x), [], pref);
    g = testclass.make(@(x) exp(-x), [], pref);
    tol_f = 10*eps*vscale(f);
    tol_g = 10*eps*vscale(g);
    pass(n, 3) = abs(innerProduct(f, g) - 2) < max(tol_f, tol_g);
    
    g = testclass.make(@(x) sin(x), [], pref);
    tol_g = 10*eps*vscale(g);
    exact = exp(1)*(sin(1) - cos(1))/2 - exp(-1)*(sin(-1) - cos(-1))/2;
    pass(n, 4) = abs(innerProduct(f, g) - exact) < max(tol_f, tol_g);
    
    %%
    % Check a few known properties.
    
    f = testclass.make(@(x) exp(x) - 1);
    g = testclass.make(@(x) 1./(1 + 1i*x.^2));
    h = testclass.make(@(x) sinh(x*exp(pi*1i/6)));
    tol_f = 10*eps*vscale(f);
    tol_g = 10*eps*vscale(g);
    tol_h = 10*eps*vscale(h);
    
    ip1 = innerProduct(alpha*f, beta*g);
    ip2 = conj(alpha)*beta*innerProduct(f, g);
    pass(n, 5) = abs(ip1 - ip2) < max(tol_f, tol_g);
    
    ip1 = innerProduct(g, h);
    ip2 = innerProduct(h, g);
    pass(n, 6) = abs(ip1 - conj(ip2)) < max(tol_g, tol_h);
    
    ip1 = innerProduct(f + g, h);
    ip2 = innerProduct(f, h) + innerProduct(g, h);
    pass(n, 7) = abs(ip1 - ip2) < max([tol_f tol_g tol_h]);
    
    ip1 = innerProduct(f, g + h);
    ip2 = innerProduct(f, g) + innerProduct(f, h);
    pass(n, 8) = abs(ip1 - ip2) < max([tol_f, tol_g, tol_h]);
    
    nf2 = innerProduct(f, f);
    ng2 = innerProduct(g, g);
    nh2 = innerProduct(h, h);
    n2vals = [nf2 ; ng2 ; nh2];
    pass(n, 9) = isreal(n2vals) && all(n2vals >= 0);
    
    %% 
    % Check operation for array-valued chebtech objects.
    
    f = testclass.make(@(x) [sin(x) cos(x)]);
    g = testclass.make(@(x) [exp(x) 1./(1 + x.^2) airy(x)]);
    tol_f = 10*max(vscale(f)*eps);
    tol_g = 10*max(vscale(g)*eps);
    ip = innerProduct(f, g);
    exact = [0.663493666631241 0                 -0.135033172317858;
             1.933421496200713 1.365866063614065  0.592109441404267];
    pass(n, 10) = norm(ip(:) - exact(:), inf) < max(tol_f, tol_g);
    
    %%
    % Check error conditions.
    
    % Can't take the inner product of a chebtech and a non-chebtech.
    try
        ip = innerProduct(f, 2); %#ok<NASGU>
        pass(n, 11) = false;
    catch ME
        pass(n, 11) = strcmp(ME.identifier, ...
            'CHEBFUN:CHEBTECH:innerProduct:input');
    end

end

end
