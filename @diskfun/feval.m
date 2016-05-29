function y = feval( f, c1, c2, c3)
%FEVAL  Evaluate a DISKFUN at one or more points.
%   Y = FEVAL( F, THETA, r)  evaluates a diskfun F at (THETA,
%   r) where THETA and r are doubles representing central angle in radians and radius on the disk.
%
%   Y = FEVAL( F, X, Y, 'cart' )  evaluates a diskfun F at a point (X,Y) in
%   Cartesian cooridnates.   
% If no coordinate flag is input, the assumption is polar coordinates. 

% See also SUBSREF.

if nargin < 4
    c3 = 1;
end

if (nargin==3) 
     theta = c1;
     r = c2; 
     y = feval@separableApprox( f, theta, r);
     if ( (size(theta, 1) == 1) && (size(r,1) == 1) )
            y = y.';
     end
elseif strcmpi(c3,'cart')  %cartesian coords are used
    if ( isnumeric(c1) && isnumeric(c2)  )
        [theta,r] = cart2pol(c1,c2); %convert to polar
        %check if any points are off the disk
        if ((any(r > 1+1e-8) ))
           error('CHEBFUN:DISKFUN:FEVAL:pointsNotOnDisk',...
            ['The specified points to evaluate the function do not '...
            'lie sufficiently close to the surface of the '...
            'unit disk.']);
        end      
        y = feval@separableApprox( f, theta, r);
        
        if ( (size(theta, 1) == 1) && (size(r,1) == 1) )
            y = y.';
        end

    elseif ( strcmp(c1, ':') && strcmp(c2, ':') )
        y = f; 
    %should we add something for evaluating slices in x and y directions?
    
    
    
    else
        error('CHEBFUN:DISKFUN:feval:argin',['Unkown input '...
            'feval(%s,%s)',c1,c2]);
    end
    
end
    

end 