function y = feval(f, c1,c2,c3)
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

%figure out if cartesian or polar

%first check global setting and apply
if strcmpi(f.coords,'cart')
    iscart =1; 
elseif strcmpi(f.coords, 'polar')
    iscart =0; 
end

%second search for user-supplied 'polar' setting in feval: 
if nargin > 3
if strcmpi(c3, 'polar')
    iscart = 0; 
end

%third search for user-supplied 'cart' setting in feval:
if ~iscart
    if strcmpi(c3, 'cart')
        iscart=1;
    end
end
end
%now evaluate
if nargin < 3 %only possibility for now is eval on a contour parametrized as complex chebfun
    c1 = imag(c1);
    c2 = real(c1); 
   % y = chebfun(@(t) feval(f, c1(t), c2(t)), 'cart') 
    y = chebfun(@(t) feval(f, real(c1(t)), imag(c1(t)), 'cart'), c1.domain, 'vectorize' );
else
    if ( isnumeric(c1) && isnumeric(c2) ) %eval at a point
        if iscart
            [theta,r] = cart2pol(c1,c2); %convert to polar
            if ((any(r > 1+1e-8) )) %check for points off disk
                error('CHEBFUN:DISKFUN:FEVAL:pointsNotOnDisk',...
                ['The specified points to evaluate the function do not '...
                'lie sufficiently close to the surface of the '...
                'unit disk.']);
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
                              c2.*sin(t)), [-pi, pi], 'trig');
    elseif (isnumeric(c1) && strcmp(c2, ':')) %radial slice
         y = chebfun(@(t) feval(f, t.*cos(c1), ...
                              t.*sin(c1)) );
    else
        error('CHEBFUN:DISKFUN:feval:argin',['Unkown input '...
            'feval(%s,%s)',c1,c2]);
    end
end

end
         
         
