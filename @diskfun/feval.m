function y = feval(varargin)
%FEVAL   Evaluate a DISKFUN at one or more points.
%   Y = FEVAL( F, X, Y) evaluates a diskfun F at a point (X,Y) in Cartesian
%   cooridnates, where X and Y are doubles.
%
%   Y = FEVAL( F, THETA, R, 'polar') evaluates a diskfun F in polar
%   coordinates (THETA,R).  Here THETA and R are doubles representing the
%   central angle (in radians) and radius in polar coordinates and must be
%   points in the unit disk.
%
%   Y = FEVAL(F, c), where c is a complex-valued chebfun representing a
%   contour, evaluates F along the contour.
%
% See also SUBSREF.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Figure out if Cartesian or polar coordinates should be used.
iscart = 1;
% Search for user-supplied 'polar' flag in arguments:
isPolar = find(strcmp(varargin, 'polar'));
if ( any( isPolar ) )
    iscart = 0;
end

% Now evaluate
f = varargin{1};
if nargin < 3 % Eval on a contour parametrized as complex chebfun
    c1 = varargin{2};
    y = chebfun(@(t) feval(f, real(c1(t)), imag(c1(t))), c1.domain, 'vectorize' );
else
    c1 = varargin{2};
    c2 = varargin{3};
    if ( isnumeric(c1) && isnumeric(c2) ) % Eval at a point
        tns =0;
        if ( ndims(c1) >= 3 && isequal(size(c1), size(c2)) )
            % x and y are tensors.
            sizec1 = size(c1);
            c1 = c1(:);
            c2 = c2(:);
            tns =1;
        end
        if iscart
            [theta,r] = cart2pol(c1,c2); % Convert to polar
            if ((any(r > 1+1e-8) )) % Check for points off disk
                error('CHEBFUN:DISKFUN:FEVAL:pointsNotOnDisk',...
                    ['The specified points to evaluate the function do not '...
                    'lie sufficiently close to the unit disk.']);
            end
            y = feval@separableApprox(f, theta, r);
        else
            theta = c1;
            r = c2;
            if ( any(r > 1+1e-8) ) % Check for points off disk
                error('CHEBFUN:DISKFUN:FEVAL:pointsNotOnDisk',...
                    ['The specified points to evaluate the function do not '...
                    'lie sufficiently close to the unit disk.']);
            end
            y = feval@separableApprox(f, theta, r);
        end

        if tns==1
            y = reshape(y, sizec1);
        end
    elseif ( strcmp(c1, ':') && strcmp(c2, ':') ) % Return the diskfun
        y = f;
    elseif ( strcmp(c1, ':') && isnumeric(c2) ) % Angular slice
        y = chebfun(@(t) feval(f, bsxfun(@times, c2, cos(t) ), ...
            bsxfun(@times, c2, sin(t) ) , 'cart'), [-pi, pi], 'trig');
    elseif (isnumeric(c1) && strcmp(c2, ':')) % Radial slice
        y = chebfun(@(t) feval(f, bsxfun(@times, cos(c1), t ),...
            bsxfun(@times, sin(c1), t ), 'cart') );
    else
        error('CHEBFUN:DISKFUN:feval:argin',['Unknown input '...
            'feval(%s,%s)',c1,c2]);
    end
end

end
