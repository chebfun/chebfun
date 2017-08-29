function g = mean2(f, dims)
%MEAN2   Average or mean value of a CHEBFUN3 in two dimensions. 
%   G = mean2(F, DIM) returns the mean of F where DIM contains two of the 
%   dimensions 1, 2 or 3 to show X, Y or Z respectively. The output is
%   a 1D CHEBFUN in the remaining variable.
%
%   G = mean2(F) is the same as mean2(F, [1,2]).
%
% See also CHEBFUN3/MEAN and CHEBFUN3/MEAN3.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check: 
if ( isempty(f) ) 
    g = []; 
    return; 
end

% Default to x and y directions: 
if ( nargin == 1 )
    dims = [1 2];
end
dim1 = dims(1); 
dim2 = dims(2);
dom = f.domain; 

if ( (dim1 == 1 && dim2 == 2) || (dim1 == 2 && dim2 == 1) )
    % mean2 over x and y: 
    dom2 = [dom(1:4)];
    
elseif ( (dim1 == 1 && dim2 == 3) || (dim1 == 3 && dim2 == 1) )
    % mean2 over x and z: 
    dom2 = [dom(1:2) dom(5:6)];
    
elseif ( (dim1 == 2 && dim2 == 3) || (dim1 == 3 && dim2 == 2) )
    % mean2 over y and z: 
    dom2 = [dom(3:6)];
    
else
    error('CHEBFUN:CHEBFUN3:mean2:dims', ...
        'dims must be [1,2], [1,3], or [2,3].');
end

g = sum2(f, dims) / (diff(dom2(1:2)) * diff(dom2(3:4)));

end