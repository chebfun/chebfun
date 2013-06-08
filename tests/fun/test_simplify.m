% Test file for fun/simplify.m

function pass = test_simplify(pref)

% Get preferences:
if ( nargin < 1 )
    pref = fun.pref;
end
pref = chebtech.pref(pref);
pref2 = pref;

pass = zeros(1, 6); % Pre-allocate pass matrix
for n = 1:1 %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();
        
        % Set the tolerance:
        tol = pref.fun.eps;
        
        % Set the domain
        dom = [-2 7];
    else
        testclass = unbndfun();
    end
    
    %%
    % Test on a scalar-valued function:
    
    f = @(x) sin(x);
    
    pref.chebtech.n = 40;
    g = testclass.make(f, dom, [], [], pref);
    h = simplify(g);
    x = h.mapping.for(h.onefun.chebpts(length(h)));
    pass(n, 1) = length(g) == 40;
    pass(n, 2) = length(h) == 25;
    pass(n, 3) = norm(f(x) - h.onefun.values, inf) < 5*h.onefun.vscale*tol;
    
    %%
    % Test on an array-valued function:
    
    f = @(x) [ sin(x), cos(x), exp(x) ];
    g = testclass.make(f, dom, [], [], pref);
    h = simplify(g);
    x = h.mapping.for(h.onefun.chebpts(length(h)));
    pass(n, 4) = length(g) == 40;
    pass(n, 5) = length(h) == 25;
    pass(n, 6) = norm(f(x) - h.onefun.values) < 5*max(h.onefun.vscale)*tol;
    
    %%
    % Test that simplifying to smaller tolerance shrinks the chebtech:
    f = testclass.make(@(x) sin(1000*(x + 0.1)), dom, [], [], pref2);
    g = simplify(f, sqrt(f.onefun.epslevel));
    pass(n, 7) = length(g) < length(f);
   
end

end
