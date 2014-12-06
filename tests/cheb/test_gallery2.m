function pass = test_gallery2(~)

names = {'airyreal', 'airycomplex', 'challenge', 'bump', 'peaks', ...
    'rosenbrock', 'smokering', 'waffle'};

N = length(names);

% Test construction of each gallery function.
for k = 1:N
    pass(k) = doesNotCrash(@() cheb.gallery2(names{k}));
end

% Test behavior for no input arguments (returning a random gallery function).
try
    pass(N+1) = doesNotCrash(@() cheb.gallery2);
end

% Test error for unknown function.
pass(N+2) = false;
try
    [f, fa] = cheb.gallery2('asdfasdfasdfasdf');
catch ME
    if strcmpi(ME.identifier, 'CHEB:GALLERY2:unknown:unknownFunction')
        pass(N+2) = true;
    end
end
end


function pass = doesNotCrash(fn)
try
    fn();
    pass = true;
catch ME %#ok<NASGU>
    pass = false;
end
end
