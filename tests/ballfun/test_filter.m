function pass = test_filter()
% Check that we can create the filter without errors

% Example 1
S = [10,10,10];
f = cheb.galleryballfun('zero',S);
pass(1) = doesNotCrash(f);

% Example 2
S = [11,12,13];
f = cheb.galleryballfun('zero',S);
pass(2) = doesNotCrash(f);

% Example 3
S = [11,13,12];
f = cheb.galleryballfun('zero',S);
pass(3) = doesNotCrash(f);

% Example 4
S = [12,11,13];
f = cheb.galleryballfun('zero',S);
pass(4) = doesNotCrash(f);

% Example 5
S = [12,13,11];
f = cheb.galleryballfun('zero',S);
pass(5) = doesNotCrash(f);

% Example 6
S = [13,11,12];
f = cheb.galleryballfun('zero',S);
pass(6) = doesNotCrash(f);

% Example 7
S = [13,12,11];
f = cheb.galleryballfun('zero',S);
pass(7) = doesNotCrash(f);

% Example 8
S = [11,11,11];
f = cheb.galleryballfun('zero',S);
pass(8) = doesNotCrash(f);

if (nargout > 0)
    pass = all(pass(:));
end
end

function pass = doesNotCrash(f)
try
    g = filter(f);  % Test returning the function
    pass = true;
catch ME %#ok<NASGU>
    pass = false;
end
end
