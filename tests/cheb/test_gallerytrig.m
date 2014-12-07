function pass = test_gallerytrig(pref)

names = {'AMsignal', 'FMsignal', 'gibbs', 'gibbsinterp', ...
    'noisyfun', 'random', 'sinefun1', 'sinefun2', 'starburst', ...
    'tsunami', 'wavepacket', 'weierstrass'};

N = length(names);

% Below, we want to test the plotting behaviour, however, we don't want the
% plots to be visible when running the test, so create an invisible figure to
% plot to:
hfig = figure('Visible', 'off');

% Test construction of each gallery function.
for k = 1:N
    pass(k) = doesNotCrash(@() cheb.gallerytrig(names{k}));
end

% Test behavior for no input arguments (returning a random gallery function).
try
    pass(N+1) = doesNotCrash(@() cheb.gallerytrig);
end

% Test error for unknown function.
pass(N+2) = false;
try
    [f, fa] = cheb.gallerytrig('asdfasdfasdfasdf');
catch ME
    if strcmpi(ME.identifier, 'CHEB:GALLERYTRIG:unknown:unknownFunction')
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
