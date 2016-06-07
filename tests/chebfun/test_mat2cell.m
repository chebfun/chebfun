% Test file for @chebfun/mat2cell.m.

function pass = test_mat2cell(pref)

% Obtain preferences.
if ( nargin == 0 )
    pref = chebfunpref();
end

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = 2 * rand(1000, 1) - 1;

% Check empty case.
C = mat2cell(chebfun());
pass(1,1) = isequal(size(C), [1 1]) && isempty(C{1});
pass(2,1) = 1;

% Check a simple example.
FF{1} = chebfun(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1], pref);
FF{2} = quasimatrix(@(x) [sin(x) cos(x) exp(x)], [-1 -0.5 0 0.5 1], pref);

% Repeat for both array-valued CHEBFUNs and quasimatrices:
for k = 1:2
    F = FF{k};

    f = chebfun(@sin, [-1 1], pref);
    g = chebfun(@cos, [-1 1], pref);
    h = chebfun(@exp, [-1 1], pref);
    fg = chebfun(@(x) [sin(x) cos(x)], pref);

    C = mat2cell(F);
    pass(k,2) = isequal(size(C), [1 3]) && (size(C{1}, 2) == 1) && ...
        (size(C{2}, 2) == 1) && (size(C{3}, 2) == 1);
    pass(k,3) = norm(feval(C{1}, xr) - feval(f, xr), inf) < ...
        10*vscale(C{1})*eps;
    pass(k,4) = norm(feval(C{2}, xr) - feval(g, xr), inf) < ...
        10*vscale(C{2})*eps;
    pass(k,5) = norm(feval(C{3}, xr) - feval(h, xr), inf) < ...
        10*vscale(C{3})*eps;

    C = mat2cell(F, [2 1]);
    pass(k,6) = isequal(size(C), [1 2]) && (size(C{1}, 2) == 2) && ...
        (size(C{2}, 2) == 1);
    err = feval(C{1}, xr) - feval(fg, xr);
    pass(k,7) = norm(err(:), inf) < 10*vscale(C{1})*eps && ...
        norm(feval(C{2}, xr) - feval(h, xr), inf) < 10*vscale(C{2})*eps;

    C = mat2cell(F, 1, [2 1]);
    pass(k,8) = isequal(size(C), [1 2]) && (size(C{1}, 2) == 2) && ...
        (size(C{2}, 2) == 1);
    err = feval(C{1}, xr) - feval(fg, xr);
    pass(k,9) = norm(err(:), inf) < 10*vscale(C{1})*eps && ...
        norm(feval(C{2}, xr) - feval(h, xr), inf) < 10*vscale(C{2})*eps;

    % Check that sizes come out right for row chebfuns.
    Ft = F.';
    ft = f.';
    gt = g.';
    ht = h.';
    fgt = fg.';

    C = mat2cell(Ft);
    pass(k,10) = isequal(size(C), [3 1]) && (size(C{1}, 1) == 1) && ...
        (size(C{2}, 1) == 1) && (size(C{3}, 1) == 1);

    C = mat2cell(Ft, [2 1]);
    pass(k,11) = isequal(size(C), [2 1]) && (size(C{1}, 1) == 2) && ...
        (size(C{2}, 1) == 1);

    C = mat2cell(Ft, [2 1], 1);
    pass(k,12) = isequal(size(C), [2 1]) && (size(C{1}, 1) == 2) && ...
        (size(C{2}, 1) == 1);

    % Check error conditions.
    try
        C = mat2cell(F, 'bad', [2 1]);
        pass(k,13) = false;
    catch ME
        pass(k,13) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:mat2cell:size');
    end

    try
        C = mat2cell(F, 2, [2 1]);
        pass(k,14) = false;
    catch ME
        pass(k,14) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:mat2cell:size');
    end

    try
        C = mat2cell(F, 1, [2 2]);
        pass(k,15) = false;
    catch ME
        pass(k,15) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:mat2cell:size');
    end

    try
        C = mat2cell(Ft, [2 1], 2);
        pass(k,16) = false;
    catch ME
        pass(k,16) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:mat2cell:size');
    end

    try
        C = mat2cell(Ft, [2 2], 1);
        pass(k,17) = false;
    catch ME
        pass(k,17) = strcmp(ME.identifier, 'CHEBFUN:CHEBFUN:mat2cell:size');
    end
end

%% Test for function defined on unbounded domain:

% Functions on [-inf b]:

% Set the domain:
dom = [-Inf -3*pi];
domCheck = [-1e6 -3*pi];

% Generate a few random points to use as test values:
x = diff(domCheck) * rand(100, 1) + domCheck(1);

% Array-valued function:
op = @(x) [exp(x) x.*exp(x) (1-exp(x))./x];
opg = @(x) exp(x);
oph = @(x) [x.*exp(x) (1-exp(x))./x];

myfun = {@chebfun, @quasimatrix};

for k = 1:2
    f = myfun{k}(op, dom);
    F = mat2cell(f, 1, [1 2]);
    F1Vals = feval(F{1}, x);
    F2Vals = feval(F{2}, x);
    
    F1Exact = opg(x);
    F2Exact = oph(x);
    err1 = F1Vals - F1Exact;
    err2 = F2Vals - F2Exact;
    
    pass(k,18) = norm([err1; err2(:)], inf) < ...
        1e2*max(eps*get(f,'vscale'));
    
end

end
