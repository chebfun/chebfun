function pass = test_gallery(~)

names = {'airy', 'bessel', 'blasius', 'bump', ...
    'chirp', 'daubechies', 'erf', 'fishfillet', 'gamma', 'gaussian', ...
    'jitter', 'kahaner', 'motto', 'random', 'rose', 'runge', ...
    'seismograph', 'Si', 'sinefun1', 'sinefun2', 'spikycomb', ...
    'stegosaurus', 'vandercheb', 'vandermonde', 'wiggly', 'wild', 'zigzag'};

N = length(names);

% Below, we want to test the plotting behaviour, however, we don't want the
% plots to be visible when running the test, so create an invisible figure to
% plot to:
hfig = figure('Visible', 'off');

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
    [f, fa] = cheb.gallery('asdfasdfasdfasdf');
catch ME
    if strcmpi(ME.identifier, 'CHEB:GALLERY:unknown:unknownFunction')
        pass(N+2) = true;
    end
end

% Close the figure used for testing the plots:
close(hfig)

end


function pass = doesNotCrash(fn)
try
    fn();           % Test plotting behavior
    [f,fa] = fn();  % Test returning the function
    pass = true;
catch ME %#ok<NASGU>
    pass = false;
end
end
