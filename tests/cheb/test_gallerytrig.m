function pass = test_gallerytrig(pref)

names = {'amsignal', 'fmsignal', 'gibbs', 'gibbsinterp', ...
    'noisyfun', 'random', 'sinefun1', 'sinefun2', 'starburst', ...
    'tsunami', 'wavepacket', 'weierstrass'};

N = length(names);

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
end


function pass = doesNotCrash(fn)
try
    fn();
    pass = true;
catch ME %#ok<NASGU>
    pass = false;
end
end
