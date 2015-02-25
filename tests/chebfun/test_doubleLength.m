function pass = test_doubleLength(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebfunpref();
end

% Test that the length doubles.
f = chebfun(@sin, pref);
g = chebfun(@sin, 'doubleLength', pref);
pass(1) = ( 2*length(f) - 1 == length(g) );

% Test on array-valued chebfuns.
f = chebfun(@(x) [sin(x) tanh(x)], pref);
g = chebfun(@(x) [sin(x) tanh(x)], 'doubleLength', pref);
pass(2) = ( 2*length(f) - 1 == length(g) );

% Test construction from values.
f = chebfun([3; 2; 1], 'doublelength');
pass(3) = (length(f) == 5);

% Test construction from coefficients.
f = chebfun([3; 2; 1], 'coeffs', 'doublelength');
pass(4) = (length(f) == 5);

% Test trig.
f = chebfun(@(x) exp(sin(pi*x)));
g = chebfun(@(x) exp(sin(pi*x)), 'doublelength');
pass(5) = (length(g) == 2*length(f) - 1);

% Test that there is an error when splitting is combined with doubleLength.
pass(6) = false;
try
    f = chebfun(@sign, 'doubleLength', 'splitting', 'on');
catch ME
    if ( strcmpi(ME.identifier, ...
            'CHEBFUN:CHEBFUN:parseInputs:doubleLengthSplitting') )
        pass(6) = true;
    end
end

% Test that there is an error when doubleLength is used with breakpoints.
pass(7) = false;
try
    f = chebfun(@exp, [0 .2 1], 'doubleLength');
catch ME
    if ( strcmpi(ME.identifier, ...
            'CHEBFUN:CHEBFUN:parseInputs:doubleLengthBreakpoints') )
        pass(7) = true;
    end
end

end
