function pass = test_doubleLength(pref)

% Get preferences:
if ( nargin < 1 )
    pref = chebfunpref();
end

% Test that the length doubles.
f = chebfun(@sin, pref);
g = chebfun(@sin, 'doubleLength', pref);
pass(1) = ( 2*length(f) + 1 == length(g) );

% Test on array-valued chebfuns.
f = chebfun(@(x) [sin(x) tanh(x)], pref);
g = chebfun(@(x) [sin(x) tanh(x)], 'doubleLength', pref);
pass(2) = ( 2*length(f) + 1 == length(g) );

% Test that there is an error when splitting is combined with doubleLength.
pass(3) = false;
try
    f = chebfun(@sign, 'doubleLength', 'splitting', 'on');
catch ME
    if ( strcmpi(ME.identifier, ...
            'CHEBFUN:CHEBFUN:parseInputs:doubleLengthSplitting') )
        pass(3) = true;
    end
end

% Test that there is an error when doubleLength is used with breakpoints.
pass(4) = false;
try
    f = chebfun(@exp, [0 .2 1], 'doubleLength');
catch ME
    if ( strcmpi(ME.identifier, ...
            'CHEBFUN:CHEBFUN:parseInputs:doubleLengthBreakpoints') )
        pass(4) = true;
    end
end

end
