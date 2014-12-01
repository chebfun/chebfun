function pass = test_gallery(pref)

names = {'airy', 'bessel', 'blasius', 'bump', 'chirp', 'erf', ...
    'fishfillet', 'gamma', 'gaussian', 'jitter', 'kahaner', 'motto', ...
    'rose', 'runge', 'seismograph', 'si', 'sinefun1', 'sinefun2', ...
    'spikycomb', 'stegosaurus', 'wiggly', 'zigzag'};

N = length(names);

% Test construction of each gallery function.
for k = 1:N
    pass(k) = doesNotCrash(@() cheb.gallery(names{k}));
end

% Test behavior for no input arguments (returning a random gallery function).
try
    pass(N+1) = doesNotCrash(@() cheb.gallery);
end

% Test error for unknown function.
pass(N+2) = false;
try
    [f,fa] = cheb.gallery('asdfasdfasdfasdf');
catch ME
    if strcmpi(ME.identifier, 'CHEB:GALLERY:unknown:unknownFunction')
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
