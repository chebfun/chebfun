function ff = mean2( f, dim )
%MEAN2   Average or mean value of a CHEBFUN3 in two dimensions.
% 
%   G = mean2(F,DIM) where DIM contains two of the dimensions 1, 2 or 3 
%   returns the mean of f only over X, Y, and Z respectively. The output is
%   a chebfun in the remaining variable.
%
%   G = mean2(F) is the same as mean2(F,[1,2])
%
%   See also MEAN and MEAN3. 

% Empty check: 
if ( isempty(f) ) 
    ff = []; 
    return; 
end

% Default to x and y directions: 
if ( nargin == 1 )
    dim = [1 2];
end
dim1 = dim(1); 
dim2 = dim(2);
dom = f.domain; 

if ( (dim1 == 1 && dim2 == 2) || (dim1 == 2 && dim2 == 1) )
    % mean2 over x and y: 
    dom2 = [dom(1:2) dom(3:4)];
    
elseif ( (dim1 == 1 && dim2 == 3) || (dim1 == 3 && dim2 == 1) )
    % mean2 over x and z: 
    dom2 = [dom(1:2) dom(5:6)];
    
elseif ( (dim1 == 2 && dim2 == 3) || (dim1 == 3 && dim2 == 2) )
    % mean2 over y and z: 
    dom2 = [dom(3:4) dom(5:6)];
    
else
    error('CHEBFUN:CHEBFUN3:sum:unknown', ...
        'Undefined function ''mean2'' for that dimension');
end

ff = sum2(f, dim) / (diff(dom2(1:2)) * diff(dom2(3:4)));

end