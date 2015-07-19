function y = feval( f, c1, c2, c3)
%FEVAL  Evaluate a DISKFUN at one or more points.
%   Y = FEVAL( F, THETA, r, 0)  evaluates a diskfun F at (THETA,
%   r) where THETA and r are doubles representing central angle in radians and radius on the disk.
%   
%   
%
%   Y = FEVAL( F, X, Y, 1 )  evaluates a diskfun F at a point (X,Y) in
%   Cartesian cooridnates.   

% See also SUBSREF.

%if nargin == 3      % Spherical coordinates used.
   % lambda = c1;
  %  theta = c2;
%elseif nargin == 4  % Cartesian coordinates used.
  %  if ~isnumeric(c1) || ~isnumeric(c1) || ~isnumeric(c1)
      %  if strcmpi(c1, ':') || strcmpi(c2, ':') || strcmpi(c3, ':') 
       %     error('SPHEREFUN:feval:colon','Colon operator not allowed when using Cartesian coordinates');
      %  else
       %     error('SPHEREFUN:feval:unknown','Unkown input feval(%s,%s,%s)',c1,c2,c3);
      %  end
    %end
    % Convert to polar coordinates
    if c3 == 1
        [theta,r] = cart2pol(c1,c2);
    else 
     theta = c1;
     r = c2; 
    % Check latitudinal coordinate system to set the elevation angle
    % appropriately.
    %if iscolat( f )
     %   theta = pi/2 - theta;
    %end
%end
    end

y = feval@separableApprox( f, theta, r);

end 