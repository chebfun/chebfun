function s = mean(f, dim)
%MEAN   Average or mean value of a CHEBFUN3.
%   MEAN(F) takes the mean of F in the x-direction (default), i.e., 
%                                           MEAN(F) = 1/(xmax-xmin) sum(F).
%
%   MEAN(F, DIM) takes the mean along the direction DIM, where DIM = 1 
%   means the x-direction, DIM = 2 means the y-direction, and DIM = 3 means 
%   the z-direction.
%
% See also CHEBFUN3/MEAN2, CHEBFUN3/MEAN3, and CHEBFUN3/STD3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(f) )
    s = chebfun3();
    return
end 

if ( nargin == 1) 
    % Default to the x-direction:
    dim = 1;  
end
dom = f.domain;

s = sum(f, dim); 
if ( dim == 1 )
    s = s / diff(dom(1:2)); % Mean in the x direction (default)
elseif ( dim == 2 )
    s = s / diff(dom(3:4)); % Mean in the y direction
elseif ( dim == 3 )
    s = s / diff(dom(5:6)); % Mean in the y direction
else
    error('CHEBFUN:CHEBFUN3:mean:dim', 'Mean not in x or y or z direction.')
end

end