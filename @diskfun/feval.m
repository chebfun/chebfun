function y = feval(varargin)
%FEVAL  Evaluate a DISKFUN at one or more points.
%   Y = FEVAL( F, THETA, r, 'polar')  evaluates a diskfun F at (THETA,
%   r) where THETA and r are doubles representing central angle in radians 
%   and radius on the disk.
%   Y = FEVAL( F, X, Y, 'cart' ) or FEVAL(F, X, Y) evaluates a diskfun F 
%   at a point (X,Y) in Cartesian cooridnates.  
%   Y = FEVAL(F, c), where c is a complex-valued chebfun representing a
%   contour, evaluates F along the contour.
%
%   See also SUBSREF.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


%figure out if cartesian or polar
 iscart = diskfun.coordsetting(varargin{:});

%now evaluate
f = varargin{1};
if nargin < 3 %eval on a contour parametrized as complex chebfun
    c1 = varargin{2};
    y = chebfun(@(t) feval(f, real(c1(t)), imag(c1(t)), 'cart'), c1.domain, 'vectorize' );
else
    c1 = varargin{2};
    c2 = varargin{3};
    if ( isnumeric(c1) && isnumeric(c2) ) %eval at a point
        if iscart
            [theta,r] = cart2pol(c1,c2); %convert to polar
            if ((any(r > 1+1e-8) )) %check for points off disk
                error('CHEBFUN:DISKFUN:FEVAL:pointsNotOnDisk',...
                ['The specified points to evaluate the function do not '...
                'lie sufficiently close to the unit disk.']);
            end 
            y = feval@separableApprox( f, theta, r);
        else
            theta = c1;
            r = c2;
        	y = feval@separableApprox( f, theta, r);
        end
        if ( (size(theta, 1) == 1) && (size(r,1) == 1) )
            y = y.';
        end
    elseif ( strcmp(c1, ':') && strcmp(c2, ':') ) %return the diskfun
        y = f;
    elseif ( strcmp(c1, ':') && isnumeric(c2) ) %angular slice
         y = chebfun(@(t) feval(f, c2.*cos(t), ...
                              c2.*sin(t), 'cart'), [-pi, pi], 'trig');
    elseif (isnumeric(c1) && strcmp(c2, ':')) %radial slice
         y = chebfun(@(t) feval(f, t.*cos(c1), ...
                              t.*sin(c1), 'cart') );
    else
        error('CHEBFUN:DISKFUN:feval:argin',['Unknown input '...
            'feval(%s,%s)',c1,c2]);
    end
end

end
         
         
