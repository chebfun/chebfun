function out = isreal( f )
%ISREAL   Real-valued CHEBFUN2 test.
%   ISREAL(F) returns logical true if F does not have an imaginary part and
%   false otherwise.
%  
%   ~ISREAL(F) detects CHEBFUN2 object that have an imaginary part even if it is
%   all zero.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty( f ) )
    out = true;
    return
end

% Get the low rank representation for f. 
cols = f.cols; 
rows = f.rows; 
piv = f.pivotValues;

% Check individual columns and rows. 
out = isreal( piv ) && isreal( cols ) && isreal( rows );

end
