function pass = test_gallery2(~)

names = {'airyreal', 'airycomplex', 'challenge', 'bump', 'peaks', ...
    'rosenbrock', 'smokering', 'waffle'};

N = length(names);

% Below, we want to test the plotting behaviour, however, we don't want the
% plots to be visible when running the test, so create an invisible figure to
% plot to:
hfig = figure('Visible', 'off');

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
