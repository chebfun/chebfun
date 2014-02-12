% Test file for fun/mat2cell.m

function pass = test_mat2cell(pref)

if ( nargin < 2 )
    pref = chebpref();
end

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
    err1 = normest(F{1} - g);
    tol1 = 10*get(g, 'epslevel')*get(g, 'vscale');
    pass(n, 1) = ~isempty(F{1}) && err1 < tol1;
    err2 = normest(F{2} - h);
    tol2 = 10*max(get(h, 'epslevel').*get(h, 'vscale'));
    pass(n, 2) = ~isempty(F{2}) && err2 < tol2;
    err3 = normest(F{3} - l);
    tol3 = 10*get(l, 'epslevel')*get(l, 'vscale');
    pass(n, 3) = ~isempty(F{3}) && err3 < tol3;
    
    % test two arguments
    F = mat2cell(f, [1 2 1]);
    err4 = normest(F{1} - g);
    tol4 = 10*get(g, 'epslevel')*get(g, 'vscale');
    pass(n, 4) = ~isempty(F{1}) && err4 < tol4;
    err5 = normest(F{2} - h);
    tol5 = 10*max(get(h, 'epslevel').*get(h, 'vscale'));
    pass(n, 5) = ~isempty(F{2}) && err5 < tol5;
    err6 = normest(F{3} - l);
    tol6 = 10*get(l, 'epslevel')*get(l, 'vscale');
    pass(n, 6) = ~isempty(F{3}) && err6 < tol6;
end

end
