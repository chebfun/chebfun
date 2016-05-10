function s = mean(f, dim)
%MEAN   Average or mean value of a SEPARABLEAPPROX. 
%   MEAN(F) takes the mean in the y-direction (default), i.e., 
%          MEAN(F) = 1/(ymax-ymin) sum(F).
%
%   MEAN(F, DIM) takes the mean along the direction DIM. If DIM = 1 it is the
%   y-direction and if DIM = 2 then it is the x-direction.
%
% See also MEAN2, STD2.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
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
    s = s / diff( dom(3:4) ); % Mean in the y direction (default)
elseif ( dim == 2 )
    s = s / diff( dom(1:2) ); % Mean in the x direction
else
    error('CHEBFUN:SEPARABLEAPPROX:mean:dim', 'Mean not in x or y direction.')
end

end
