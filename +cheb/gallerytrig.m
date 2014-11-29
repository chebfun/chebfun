function [f,fa] = gallerytrig(name)
%GALLERYTRIG   Chebfun periodic test functions.
%   GALLERYTRIG(NAME) returns a periodic chebfun corresponding to NAME. See
%   the listing below for available function names.
%
%   [F,FA] = GALLERYTRIG(NAME) also returns the anonymous function used to
%   define the function.
%
%   gibbs        Gibbs-Wilbraham approximation of a square wave
%   gibbsinterp  Interpolant of a square wave
%   weierstrass  The first eight terms of the pathological function

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


% If the user did not supply an input, return a random function
% from the gallery.
if ( nargin == 0 )
    names = {'gibbs', 'gibbsinterp', 'weierstrass'};
    name = names{randi(length(names))};
end

% The main switch statement.
switch name

    % Function from SIAM 100-digit challenge
    case 'gibbs'
        n = 31;
        N = (1:2:n)';
        fa = @(x) sum(4/pi * sin(bsxfun(@times,N,x(:)')) ...
                            ./ (N*ones(1,length(x))) )';
        f = chebfun(fa, [-pi pi], 'trig');

    % Interpolant of a square wave
    case 'gibbsinterp'
        n = 31;
        fa = @(x) 2*(x > 0) - 1;
        f = chebfun(fa, [-pi pi], 'trig', 2*n+1);

    % The first eight terms of the pathological function
    case 'weierstrass'
        k = 8;
        K = (1:k)';
        fa = @(x) sum( 2.^-(K*ones(1,length(x))) ...
                    .* cos(bsxfun(@times, 4.^K, x(:)')) )';
        f = chebfun(fa, [-pi/4 pi/4], 'trig');

    % Error if the input is unknown.
    otherwise
        error('CHEB:GALLERY:unknown:unknownFunction', ...
            'Unknown function.')

end

end
