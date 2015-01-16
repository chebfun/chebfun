function pass = test_bspline(~)
%TEST_BERNOULLI   Test that the cheb.bpsline method doesn't crash.

pass = doesNotCrash(@() cheb.bspline(4));

end

function pass = doesNotCrash(fn)
try
    fn();
    pass = true;
catch ME %#ok<NASGU>
    pass = false;
end
end
