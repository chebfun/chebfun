function f = ctranspose(f)
%'	  Complex conjugate transpose of a chebfun2.
% 
% F' is the complex conjugate transpose of F.
%
% G = CTRANSPOSE(F) is called for the syntax F'.  

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Empty check
if ( isempty( f ) )
    f = chebfun2; 
    return
end

% Get the low rank representation for f. 
cols = f.cols; 
rows = f.rows; 

% Swap columns and rows:
temp = cols; 
f.cols = rows; 
f.rows = temp; 

end