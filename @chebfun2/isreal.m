function isReal=isreal( f )
%ISREAL Real-valued chebfun2 test.
%
% ISREAL(F) returns logical true if F does not have an imaginary part
% and false otherwise.
%  
% ~ISREAL(F) detects chebfun2s that have an imaginary part even if
% it is all zero.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

if ( isempty( f ) )
    isReal = [];
    return
end

% Get the low rank representation for f. 
cols = f.cols; 
rows = f.rows; 
piv = f.pivotValues;

% Check individual columns and rows. 
isReal = isreal( cols ) && isreal( rows ) && isreal( piv );

end
