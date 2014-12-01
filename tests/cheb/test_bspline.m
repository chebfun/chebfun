function pass = test_bspline(pref)

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
