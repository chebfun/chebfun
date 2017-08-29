function varargout = gallery2(name)
%CHEB.GALLERY2   Chebfun2 example functions.
%   F = CHEB.GALLERY2(NAME) returns a chebfun2 corresponding to NAME. See the
%   listing below for available names.
%
%   For example,  plot(cheb.gallery2('peaks'))  plots the classic MATLAB peaks
%   function represented as a chebfun2. For details of how each function is
%   constructed, try type +cheb/gallery2 or edit cheb.gallery2.
%
%   [F,FA] = CHEB.GALLERY2(NAME) also returns the anonymous function FA used to
%   define the function.
%
%   CHEB.GALLERY2 with no input argument returns a function selected at
%   random from the gallery.
%
%   CHEB.GALLERY2 with no output argument creates a plot of the selected
%   function.
%
%   airyreal     Real part of Airy Ai function
%   airycomplex  Imaginary part of Airy Ai function
%   bump         2D C-infinity function with compact support
%   challenge    Function from SIAM 100-digit challenge
%   peaks        Classic MATLAB peaks function
%   rosenbrock   Challenging test function in optimization
%   roundpeg     Approx characteristic function of a disk, rank 45
%   smokering    A halo, hoop, or hole
%   squarepeg    Approx characteristic function of a square, rank 1
%   tiltedpeg    A tilted variant of squarepeg, rank 100
%   waffle       Function with horizontal and vertical ridges
%
%   Gallery functions are subject to change in future releases of Chebfun.
%
% See also CHEB.GALLERY, CHEB.GALLERYTRIG, CHEB.GALLERY3, CHEB.GALLERYDISK, CHEB.GALLERYSPHERE.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% If the user did not supply an input, return a function selected at
% random from the gallery.
if ( nargin == 0 )
    names = {'airyreal', 'airycomplex', 'challenge', 'bump', 'peaks', ...
        'rosenbrock', 'roundpeg', 'smokering', 'squarepeg', 'tiltedpeg', ...
        'waffle'};
    name = names{randi(length(names))};
end

% The main switch statement.
switch lower(name)

    % Real part of Airy Ai function:
    case 'airyreal'
        fa = @(z) real(airy(z));
        f = chebfun2(fa, [-10 10 0 1]);

    % Imaginary part of Airy Ai function:
    case 'airycomplex'
        fa = @(z) airy(z);
        f = chebfun2(fa, [-10 10 -5 5]);

    % Two-dimensional bump function:
    case 'bump'
        fa = @(x,y) (x.^2 + y.^2 < 1).*exp(-(x.^2 + y.^2 < 1)./ ...
            (1 - x.^2 - y.^2));
        f = chebfun2(fa, [-1 1 -1 1]*2);

    % Function from SIAM 100-digit challenge:
    case 'challenge'
        fa = @(x,y) exp(sin(50*x)) + sin(60*exp(y)) + sin(70*sin(x)) + ...
            sin(sin(80*y)) - sin(10*(x+y)) + (x.^2+y.^2)./4;
        f = chebfun2(fa);

    % The classic MATLAB peaks function:
    case 'peaks'
        fa = @(x,y) 3*(1-x).^2.*exp(-x.^2 - (y+1).^2) ...
            - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2 - y.^2) ...
            - 1/3*exp(-(x+1).^2 - y.^2);
        f = chebfun2(fa, [-3 3 -3 3]);

    % A challenging test function in optimization:
    case 'rosenbrock'
        fa = @(x,y) (1-x).^2 + 100*(y-x.^2).^2;
        f = chebfun2(fa, [-2 2 -1 3]);

    % Approximately the characteristic function of a disk:
    case 'roundpeg'
        fa = @(x,y) 1./(1+((2*x).^2+(2*y).^2).^10);
        f = chebfun2(fa);
        
    % A halo, hoop, or hole:
    case 'smokering'
        fa = @(x,y) exp(-100*(x.^2 - x.*y + 2*y.^2 - 1/2).^2);
        f = chebfun2(fa);

    % Approximately the characteristic function of a square:
    case 'squarepeg'
        fa = @(x,y) 1./((1+(2*x).^20).*(1+(2*y).^20));
        f = chebfun2(fa);

    % A tilted version of squarepeg:
    case 'tiltedpeg'
        fa = @(x,y) 1./((1+(2*x+.4*y).^20).*(1+(2*y-.4*x).^20));
        f = chebfun2(fa);

    % A function with horizontal and vertical ridges:
    case 'waffle'
        fa = @(x,y) 1./(1 + 1e3*((x.^2-.25).^2.*(y.^2-.25).^2));
        f = chebfun2(fa);

    % Raise an error if the input is unknown.
    otherwise
        error('CHEB:GALLERY2:unknown:unknownFunction', ...
            'Unknown function.')
end

% Only return something if there is an output argument.
if ( nargout > 0 )
    varargout = {f, fa};
else
    % Otherwise, plot the function.
    plot(f)
    title([name ', rank = ' num2str(length(f))])
end

end
