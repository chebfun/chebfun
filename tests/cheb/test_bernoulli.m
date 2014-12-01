function pass = test_bernoulli(~)
%TEST_BERNOULLI   Test that the Bernoulli method doesn't crash.
pass = doesNotCrash(@() cheb.bernoulli(4));

end

function pass = doesNotCrash(fn)
try
    fn();
    pass = true;
catch ME %#ok<NASGU>
    pass = false;
end
end
