function s = mean2(f, dims)
%MEAN   Average or mean value of a BALLFUN in two directions. 
%   MEAN(F, DIM) computes the mean of F over two of the variables r, lambda or theta 
%   where DIMS is a row vector containing two of the three indices
%   1, 2 or 3. The output is a 1D CHEBFUN in the remaining variable.
%
% See also MEAN2, MEAN3. 

% Copyright 2019 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Default to lambda and theta directions: 
if ( nargin == 1 )
    dims = [2 3];
end

s = sum2(f, dims);

% Sort the dimensions
dims = sort(dims);
dim_1 = dims(1); dim_2 = dims(2);

if dim_1 == 2 && dim_2 == 3
    % Average over lambda and theta
    s = s/(4*pi);
    
elseif dim_1 == 1 &&  dim_2 == 2
    % Average over r and lambda
    s = 3*s/(2*pi);

elseif dim_1 == 1 && dim_2 == 3
    % Average over r and theta
    s = 3*s/2;
else
    error('CHEBFUN:BALLFUN:mean2:dim', ...
        'Unrecognized input.')
end    
end