% Test file for fun/mat2cell.m

function pass = test_mat2cell(pref)

if ( nargin < 2 )
    pref = fun.pref;
end

pass = zeros(1, 2); % Pre-allocate pass matrix.
for n = 1:1 %[TODO]: unbndfun
    if ( n == 1 )
        testclass = bndfun();
        dom = [-2 7];
    else 
        testclass = unbndfun();
    end

    f = testclass.make(@(x) [sin(x) cos(x) exp(x) x], dom, [], [], pref);
    g = testclass.make(@(x) sin(x), dom, [], [], pref);
    h = testclass.make(@(x) [cos(x) exp(x)], dom, [], [], pref);
    l = testclass.make(@(x) x, dom, [], [], pref);
    
    % test full arguments
    F = mat2cell(f, 1, [1 2 1]);
    pass(n, 1) = sum(F{1} - g) < g.onefun.epslevel;
    pass(n, 2) = all( sum(F{2} - h) < h.onefun.vscale*h.onefun.epslevel );
    pass(n, 3) = sum(F{3} - l) < l.onefun.epslevel;
    
    % test two arguments
    F = mat2cell(f, [1 2 1]);
    pass(n, 4) = sum(F{1} - g) < g.onefun.epslevel;
    pass(n, 5) = all( sum(F{2} - h) < h.onefun.vscale*h.onefun.epslevel );
    pass(n, 6) = sum(F{3} - l) < l.onefun.epslevel;
end

end
