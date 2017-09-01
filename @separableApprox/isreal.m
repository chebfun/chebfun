function out = isreal( f )
%ISREAL   Real-valued SEPARABLEAPPROX test.
%   ISREAL(F) returns logical true if F does not have an imaginary part and
%   false otherwise.
%  
%   ~ISREAL(F) detects SEPARABLEAPPROX object that have an imaginary part even if it is
%   all zero.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

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
