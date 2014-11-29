function [f,fa] = gallery(name)
%GALLERY   Chebfun test functions.
%   GALLERY(NAME) returns a chebfun corresponding to NAME. See the listing
%   below for available function names.
%
%   [F,FA] = GALLERY(NAME) also returns the anonymous function used to define
%   the function.
%
%   airy         Airy Ai function on [-40, 40]
%   bessel       Bessel function with parameter 0 on [-100, 100]
%   bump         A bump function
%   blasius      The Blasius function on [0, 10]
%   chirp        Sine with exponentially increasing frequency
%   erf          The error function on [-10, 10]
%   fishfilet    Wild oscillations from Extreme Extrema example
%   gamma        The gamma function on [-4, 4]
%   jitter       A piecewise constant function
%   kahaner      Challenging integrand with four spikes
%   motto        (Scribbled) Chebfun motto by Gilbert Strang
%   runge        The Runge function
%   seismograph  Tanh plus growing oscillation
%   si           Sine integral on [-50, 50]
%   sinefun1     As smooth as it looks
%   sinefun2     Not as smooth as it looks
%   spikycomb    25 peaks, each sharper than the last
%   wiggly       One of the Chebfun team's favorites
%   zigzag       Degree 10000 polynomial that looks piecewise linear

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% If the user did not supply an input, return a random function
% from the gallery.
if ( nargin == 0 )
    names = {'airy', 'bessel', 'chirp', 'erf', 'fishfilet', 'gamma', ...
            'jitter', 'kahaner', 'motto', 'runge', 'sinefun1', 'sinefun2', ...
            'si', 'spikycomb', 'seismograph', 'wiggly', 'zigzag'};
    name = names{randi(length(names))};
end

% The main switch statement.
switch name

    % The Airy function on [-40, 40]
    case 'airy'
        fa = @airy;
        f = chebfun(fa, [-40 40]);

    % Bessel function with parameter 0 on [-100, 100]
    case 'bessel'
        fa = @(x) besselj(0, x);
        f = chebfun(fa, [-100 100]);

    % The Blasius function on [0, 10]
    case 'blasius'
        op = @(u) 2*diff(u,3) + u.*diff(u,2);
        bc  = @(x,u) [u(0); feval(diff(u),0); feval(diff(u),10)-1];
        N = chebop(op, [0 10], bc);
        N.init = chebfun([3.55542; 4.25907; 0.43669; -0.21367;
                0.06382; -0.00118; -0.00865; 0.00306], [0 10], 'coeffs');
        f = N\0;
        fa = @(x) f(x);

    % A bump function
    case 'bump'
        fa = @(x) (abs(x) < 1).*exp(-1./(1-x.^2));
        f = chebfun(fa, [-2 2]);

    % Sine with exponentially increasing frequency
    case 'chirp'
        fa = @(x) sin(x.*exp(x));
        f = chebfun(fa, [0 5]);

    % The error function on [-10, 10]
    case 'erf'
        fa = @erf;
        f = chebfun(fa, [-10 10]);

    % Wild oscillations from Extreme Extrema example
    case 'fishfilet'
        fa = @(x) cos(x).*sin(exp(x));
        f = chebfun(fa, [0 6]);

    % The gamma function on [-4, 4]
    case 'gamma'
        fa = @gamma;
        f = chebfun(fa, [-4 4], 'blowup', 'on', 'splitting', 'on');

    % A piecewise constant function
    case 'jitter'
        fa = @(x) round(exp(x)*2.*sin(8*x));
        f = chebfun(fa, 'splitting', 'on');

    % Challenging integrand with four spikes
    case 'kahaner'
        fa = @(x) sech(10*(x-0.2)).^2 + sech(100*(x-0.4)).^4 + ...
         sech(1000*(x-0.6)).^6 + sech(1000*(x-0.8)).^8;
        f = chebfun(fa, [0 1]);

    % (Scribbled) Chebfun motto by Gilbert Strang
    case 'motto'
        f = exp(3i*scribble('there is no fun like chebfun'));
        fa = @(x) f(x);

    case 'rose'
        m = 5; n = 4;
        fa = @(t) cos(m/n*t).*cos(t) + 1i*cos(m/n*t).*sin(t);
        f = chebfun(fa, [0, 8*pi], 'trig');

    % The Runge function
    case 'runge'
        fa = @(x) 1./(1 + x.^2);
        f = chebfun(fa, [-5 5]);

    % As smooth as it looks
    case 'sinefun1'
        fa = @(x) (2 + sin(50*x));
        f = chebfun(fa);

    % Not as smooth as it looks
    case 'sinefun2'
        fa = @(x) (2 + sin(50*x)).^1.0001;
        f = chebfun(fa);

    % Sine integral
    case 'si'
        f = cumsum(chebfun(@(x) sin(x)./(x)), [-50, 50]);
        fa = @(x) f(x);

    % 25 peaks, each sharper than the last
    % from ATAP, Chapter 18
    case 'spikycomb'
        fa = @(x) exp(x).*sech(4*sin(40*x)).^exp(x);
        f = chebfun(fa);

    % Tanh plus growing oscillation
    % from ATAP, Chapter 5
    case 'seismograph'
        fa = @(x) tanh(20*sin(12*x)) + .02*exp(3*x).*sin(300*x);
        f = chebfun(fa);

    % One of the Chebfun team's favorites
    case 'wiggly'
        fa = @(x) sin(x) + sin(x.^2);
        f = chebfun(fa, [0 10]);

    % Degree 10000 polynomial that looks piecewise linear
    % from ATAP appendix
    case 'zigzag'
        f = cumsum(chebfun(@(t) sign(sin(100*t./(2-t))), 10000));
        fa = @(x) f(x);

    % Error if the input is unknown.
    otherwise
        error('CHEB:GALLERY:unknown:unknownFunction', ...
            'Unknown function.')

end

end
