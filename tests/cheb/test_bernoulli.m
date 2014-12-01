function pass = test_bernoulli(pref)

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
