% Test file for chebtech/trigcoeffs.m

function pass = test_trigcoeffs(pref)

% Get preferences.
if ( nargin < 1 )
    pref = chebtech.techPref();
end

for n = 1:2
    if ( n == 1 )
        testclass = chebtech1();
    else 
        testclass = chebtech2();
    end
    
    %%
    % Check a few simple examples.

    f = testclass.make(@(x) zeros(size(x)), [], pref);
    p = trigcoeffs(f);
    pass(n,1) = (norm(p, inf) <= 10*f.vscale.*f.epslevel);

    f = testclass.make(@(x) 3*ones(size(x)), [], pref);
    p = trigcoeffs(f);
    pass(n,2) = (norm(p - 3, inf) < 10*f.vscale.*f.epslevel);

    f = testclass.make(@(x) 1+cos(pi*x), [], pref);
    p = trigcoeffs(f,3);
    pass(n,3) = (norm(p - [0.5 1 0.5]', inf) < 10*f.vscale.*f.epslevel);
    p = trigcoeffs(f,5);
    pass(n,4) = (norm(p - [0 0.5 1 0.5 0]', inf) < 10*f.vscale.*f.epslevel);
    p = trigcoeffs(f,1);
    pass(n,5) = (norm(p - 1, inf) < 10*f.vscale.*f.epslevel);

    f = testclass.make(@(x) 1 + exp(2*1i*pi*x) + exp(-1i*pi*x), [], pref);
    p = trigcoeffs(f,5);
    pass(n,6) = (norm(p - [0 1 1 0 1 ]', inf) ...
        < 10*f.vscale.*f.epslevel);
    p = trigcoeffs(f,9);
    pass(n,7) = (norm(p - [0 0 0 1 1 0 1 0 0 ]', inf) ...
        < 10*f.vscale.*f.epslevel);
    p = trigcoeffs(f,3);
    pass(n,8) = (norm(p - [1 1 0]', inf) ...
        < 10*f.vscale.*f.epslevel);

    %%
    % Verify operation for array-valued chebtech objects.

    f = testclass.make(@(x) [3*ones(size(x)), 1+cos(pi*x), ... 
        1 + exp(2*1i*pi*x) + exp(-1i*pi*x)], [], pref);
    p = trigcoeffs(f,5);
    p_exact = [0 0   0;...
               0 0.5 1;...   
               3 1   1;...
               0 0.5 0;...
               0 0   1];
    pass(n,9) = (norm(p(:) - p_exact(:), inf) < 10*max(f.vscale.*f.epslevel));

    p = trigcoeffs(f,7);
    p_exact = [0 0   0;...
               0 0   0;...
               0 0.5 1;...   
               3 1   1;...
               0 0.5 0;...
               0 0   1;...
               0 0   0];
    pass(n,10) = (norm(p(:) - p_exact(:), inf) < 10*max(f.vscale.*f.epslevel));

    p = trigcoeffs(f,3);
    p_exact = [0 0.5 1;...   
               3 1   1;...
               0 0.5 0];
    pass(n,11) = (norm(p(:) - p_exact(:), inf) < 10*max(f.vscale.*f.epslevel));

    p = trigcoeffs(f,0);
    pass(n,12) = isempty(p);
end

end
