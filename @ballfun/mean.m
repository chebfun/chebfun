function s = mean(f, dim)
%MEAN   Average or mean value of a BALLFUN in a specific direction. 
%   MEAN(F, DIM) where DIM is 1, 2 or 3 is the mean of F over r (radial direction), 
%   lambda (azimuthal direction) or theta (polar direction) respectively and 
%   and returns as its output a spherefun if DIM is 1 or a diskfun otherwise.
%
% See also BALLFUN/MEAN2, BALLFUN/MEAN3. 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( nargin == 1) 
    % Default to the r-direction:
    dim = 1;    
end

% Compute the definite integral of f
s = sum(f, dim);

if dim == 1
    % Mean in the r direction
    s = s*3;
elseif dim == 2
    % Mean in the lambda direction
    s = s/(2*pi);
elseif dim == 3
    % Mean in the theta direction
    s = s/2;
else
    error('CHEBFUN:BALLFUN:mean:dim', ...
        'Unrecognized input.')
end    
end
