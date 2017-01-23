function s = mean(f, dim)
%MEAN   Average or mean value of a SPHEREFUN. 
%   MEAN(F) takes the mean in the latitude-direction (default), i.e., 
%          MEAN(F) = 1/pi sum(F).
%
%   MEAN(F, DIM) takes the mean along the direction DIM. If DIM = 1 it is the
%   latitude-direction and if DIM = 2 then it is the longitude-direction.
%
% See also SPHEREFUN/MEAN2, SPHEREFUN/STD2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( f ) )
    s = chebfun;
    return
end 

if ( nargin == 1) 
    % Default to the y-direction:
    dim = 1;    
end
dom = f.domain;

s = sum( f, dim ); 
if ( dim == 1 )
    s = s / pi; % Mean in the latitude direction (default)
elseif ( dim == 2 )
    s = s / (2*pi); % Mean in the longitude direction
else
    error('CHEBFUN:SPHEREFUN:mean:dim', ...
        'Mean not in longitude (LAMBDA) or latitude (THETA) direction.')
end

end

