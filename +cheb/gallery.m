function [f,fa] = gallery(name)
%GALLERY   Gallery of 1-dimensional functions.
%   GALLERY(N) returns interesting 1D functions as a CHEBFUN quasimatrix.
%   N must be a vector with integer entries taking values from 1 to 11.
%   All gallery functions have domain [-1, 1].
%
%   [F,FA] = GALLERY(N) also returns the anonymous functions used to define
%   the CHEBFUN. If N is an integer, FA is a function handle; if N is a
%   vector, FA is a cell array of function handles.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

switch name

    % the Airy function on [-40 40]
    case 'airy'
        fa = @airy;
        f = chebfun(fa, [-40 40]);

    % the airy function on [-40 40]
    case {'bessel', 'besselj'}
        fa = @(x) besselj(0, x);
        f = chebfun(fa, [-100 100]);

    % sine with exponentially increasing frequency
    case 'chirp'
        fa = @(x) sin(x.*exp(x));
        f = chebfun(fa, [0 5]);

    % the error function
    case 'erf'
        fa = @erf;
        f = chebfun(fa, [-10 10]);

    % wild oscillations from Extreme Extrema example
    case 'fishfilet'
        fa = @(x) cos(x).*sin(exp(x));
        f = chebfun(fa, [0 6]);

    % the gamma function on [-4, 4]
    case 'gamma'
        fa = @gamma;
        f = chebfun(fa, [-4 4], 'blowup', 'on', 'splitting', 'on');

    % a piecewise constant function
    case 'jitter'
        fa = @(x) round(exp(x)*2.*sin(8*x));
        f = chebfun(fa, 'splitting', 'on');

    % challenging integrand with four spikes
    case 'kahaner'
        fa = @(x) sech(10*(x-0.2)).^2 + sech(100*(x-0.4)).^4 + ...
         sech(1000*(x-0.6)).^6 + sech(1000*(x-0.8)).^8;
        f = chebfun(fa, [0 1]);

    % (scribbled) Chebfun motto by Gilbert Strang
    case 'motto'
        f = exp(3i*scribble('there is no fun like chebfun'));
        fa = @(x) f(x);

    % The Runge function
    case 'runge'
        fa = @(x) 1./(1 + x.^2);
        f = chebfun(fa, [-5 5]);

    % as smooth as it looks
    case 'sinefun1'
        fa = @(x) (2 + sin(50*x));
        f = chebfun(fa);

    % not as smooth as it looks
    case 'sinefun2'
        fa = @(x) (2 + sin(50*x)).^1.0001;
        f = chebfun(fa);

    % degree 10000 polynomial that looks piecewise linear
    case 'si'
        f = cumsum(chebfun(@(x) sin(50*x)./(50*x)));
        fa = @(x) f(x);

    % 25 peaks, each sharper than the last
    % from ATAP, Chapter 18
    case 'spikycomb'
        fa = @(x) exp(x).*sech(4*sin(40*x)).^exp(x);
        f = chebfun(fa);

    % tanh plus growing oscillation
    % from ATAP, Chapter 5
    case 'wiggles'
        fa = @(x) tanh(20*sin(12*x)) + .02*exp(3*x).*sin(300*x);
        f = chebfun(fa);

    % one of the Chebfun team's favorites
    case 'wiggly'
        fa = @(x) sin(x) + sin(x.^2);
        f = chebfun(fa, [0 10]);

    % degree 10000 polynomial that looks piecewise linear
    % from ATAP appendix
    case 'zigzag'
        f = cumsum(chebfun(@(t) sign(sin(100*t./(2-t))), 10000));
        fa = @(x) f(x);

end

end
