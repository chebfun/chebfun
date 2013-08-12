function pass = test_constructor_inputs(pref)

if ( nargin == 0 )
    pref = chebfun.pref();
end

% [TODO]: This test needs to be updated to include more exotic input options.

% Test the vectorise flag:
try
    f = chebfun(@(x) x^2, pref, 'vectorize'); %#ok<NASGU>
    pass(1) = true;
catch
    pass(1) = false;
end

% Test the vectorise flag:
f = chebfun(@(x) sin(x), pref, 5);
pass(2) = ( length(f.funs{1}) == 5 ); 
% [TODO]: Change this once @CHEBFUN/LENGTH is implemented.

% Test the 'splitting on' flag.
f = chebfun(@(x) abs(x), 'splitting', 'on');
pass(3) = numel(f.funs) == 2;

% Test the 'extrapolate' flag.
f = chebfun(@(x) sign(x), [0, 1], 'extrapolate', 'on');
pass(4) = get(f, 'ishappy');

end