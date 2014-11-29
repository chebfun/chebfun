function [f,fa] = gallery2(name)
%GALLERY2   Chebfun2 test functions.
%   GALLERY2(NAME) returns a chebfun2 corresponding to NAME. See the listing
%   below for available function names.
%
%   [F,FA] = GALLERY2(NAME) also returns the anonymous function used to define
%   the test function.
%
%   bump         A two-dimensional bump function
%   challenge    Function from SIAM 100-digit challenge
%   smokering    A halo, hoop, or hole.
%   peaks        The classic MATLAB peaks function
%   rosenbrock   A challenging test function in optimization
%   waffle       A function with horizontal and vertical ridges

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% If the user did not supply an input, return a random function
% from the gallery.
if ( nargin == 0 )
    names = {'challenge', 'bump', 'peaks', 'rosenbrock', 'smokering', ...
            'waffle'};
    name = names{randi(length(names))};
end

% The main switch statement.
switch name

    % Function from SIAM 100-digit challenge
    case 'challenge'
        fa = @(x,y) exp(sin(50*x)) + sin(60*exp(y)) + sin(70*sin(x)) + ...
            sin(sin(80*y)) - sin(10*(x+y)) + (x.^2+y.^2)./4;
        f = chebfun2(fa);

    % A two-dimensional bump function
    case 'bump'
        fa = @(x,y) (x.^2 + y.^2 < 1).*exp(-(x.^2 + y.^2 < 1)./(1 - x.^2 - y.^2));
        f = chebfun2(fa, [-1 1 -1 1]*2);

    % The classic MATLAB peaks function
    case 'peaks'
        fa = @(x,y) 3*(1-x).^2.*exp(-x.^2 - (y+1).^2) ...
            - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2 - y.^2) ...
            - 1/3*exp(-(x+1).^2 - y.^2);
        f = chebfun2(fa, [-3 3 -3 3]);

    % A challenging test function in optimization
    case 'rosenbrock'
        fa = @(x,y) (1-x).^2 + 100*(y-x.^2).^2;
        f = chebfun2(fa, [-2 2 -1 3]);

    % A halo, hoop, or hole.
    case 'smokering'
        fa = @(x,y) exp(-100*(x.^2-x.*y+2*y.^2 - 1/2).^2);
        f = chebfun2(fa);

    % A function with horizontal and vertical ridges
    case 'waffle'
        fa = @(x,y) 1./(1+1e3*((x.^2-.25).^2.*(y.^2-.25).^2));
        f = chebfun2(fa);

    % Error if the input is unknown.
    otherwise
        error('CHEB:GALLERY:unknown:unknownFunction', ...
            'Unknown function.')

end

end
