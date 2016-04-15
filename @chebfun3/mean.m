function s = mean(f, dim)
%MEAN   Average or mean value of a CHEBFUN3. 
%   MEAN(F) takes the mean in the x-direction (default), i.e., 
%          MEAN(F) = 1/(xmax-xmin) sum(F).
%
%   MEAN(F, DIM) takes the mean along the direction DIM. If DIM = 1 it is the
%   x-direction, and if DIM = 2 then it is the y-direction, and if DIM = 3 then it is the z-direction.
%
%   See also MEAN2, MEAN3, and STD3.

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