% Test file for @chebtech/sample.m

function pass = test_sample(pref)

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

    f = testclass.make(@(x) sin(x - 0.1));

    % Test on a grid equal to length(f).
    [v, p] = sample(f);
    p_ex = testclass.chebpts(length(f));
    v_ex = feval(f, p_ex);
    pass(n, 1) = (norm(p - p_ex) < 100*eps) && (norm(v - v_ex) < 100*eps);

    % Test on a grid shorter than length(f).
    m = round(length(f)/2);
    [v, p] = sample(f, m);
    p_ex = testclass.chebpts(m);
    v_ex = feval(f, p_ex);
    pass(n, 2) = (norm(p - p_ex) < 100*eps) && (norm(v - v_ex) < 100*eps);

    % Test on a grid longer than length(f).
    m = round(2*length(f));
    [v, p] = sample(f, m);
    p_ex = testclass.chebpts(m);
    v_ex = feval(f, p_ex);
    pass(n, 3) = (norm(p - p_ex) < 100*eps) && (norm(v - v_ex) < 100*eps);
end

end
