function pass = test_isnan(pref)

if ( nargin < 1 )
    pref = funcheb.pref;
end
p = pref;

% Test a scalar-valued function:
f = funcheb2(@(x) x, p);
pass(1) = ~isnan(f);

% Test a vector-valued function:
f = funcheb2(@(x) [x, x.^2], p);
pass(2) = ~isnan(f);

% Test a NaN scalar-valued function:
try
    f = funcheb2(@(x) x + NaN);
    pass(3) = isnan(f);
catch ME
    pass(3) = strcmpi(ME.message, 'Too many NaNs to handle.');
end

% Test a NaN vector-valued function:
try
    f = funcheb2(@(x) [x, x + NaN]);
    pass(4) = isnan(f);
catch ME
    pass(4) = strcmpi(ME.message, 'Too many NaNs to handle.');
end

% Test a NaN vector-valued function:
try
    f = funcheb2(@(x) myfun(x))
    pass(5) = isnan(f);
catch ME
    pass(5) = strcmpi(ME.message, 'Function returned NaN when evaluated.');
end

% Test a non-adaptive construction
p.funcheb.n = 11;
try
    f = funcheb2(@(x) myfun(x), p);
    pass(6) = isnan(f);
catch ME
    pass(6) = strcmpi(ME.message, 'Function returned NaN when evaluated.');
end

end


function y = myfun(x)
    y = 1./x;
    y(x == 0) = NaN;
end