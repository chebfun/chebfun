function f = transpose( f ) 
% .'   SEPARABLEAPPROX transpose. 
%    F.' is the non-conjugate transpose of a F. 
% 
% See also CTRANSPOSE. 

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Empty check:
if ( isempty( f ) )
    return
end

% Swap columns and rows. 
temp = f.cols; 
f.cols = f.rows; 
f.rows = temp; 

end
