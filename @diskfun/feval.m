function y = feval( f, c1, c2, c3)
%FEVAL  Evaluate a DISKFUN at one or more points.
%   Y = FEVAL( F, THETA, r, 1)  evaluates a diskfun F at (THETA,
%   r) where THETA and r are doubles representing central angle in radians and radius on the disk.
%
%   Y = FEVAL( F, X, Y, 0 )  evaluates a diskfun F at a point (X,Y) in
%   Cartesian cooridnates.   
% If no coordinate flag is input, the assumption is polar coordinates. 

% See also SUBSREF.

if nargin < 4
    c3 = 1;
end

    % Convert to polar coordinates
    if c3 == 0
        [theta,r] = cart2pol(c1,c2);
    else 
     theta = c1;
     r = c2; 
    end

y = feval@separableApprox( f, theta, r);

end 