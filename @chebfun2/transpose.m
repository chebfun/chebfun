function f = transpose( f ) 
% .'   CHEBFUN2 transpose. 
%    F.' is the non-conjugate transpose of a F. 
% 
% See also CTRANSPOSE. 

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check:
if ( isempty( f ) )
    return
end

% Swap columns and rows. 
temp = f.cols; 
f.cols = f.rows; 
f.rows = temp; 

end
