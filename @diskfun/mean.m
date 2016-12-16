function s = mean(f, dim)
%MEAN   Average or mean value of a DISKFUN. 
%   MEAN(F) takes the mean in the angular direction (default), i.e., 
%          MEAN(F) = sum(F).
%
%   MEAN(F, DIM) takes the mean along the direction DIM. If DIM = 1 it is the
%   radial direction, and if DIM = 2 then it is the angular direction,
%   i.e., MEAN(F,2) = 1/(2*pi)*sum(F,2).
%
% See also DISKFUN/MEAN2, DISKFUN/STD2.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty(f) )
    s = chebfun;
    return
end 

if ( nargin == 1)
    % Default to the angular direction.
    dim = 1;    
end

s = sum(f, dim);
if ( dim == 1 )
    return; % Mean in the angular direction (default)
elseif ( dim == 2 )
    s = s / (2*pi); % Mean in the radial direction
else
    error('CHEBFUN:DISKFUN:mean:dim', ...
        'dim must be 1 or 2.')
end

end