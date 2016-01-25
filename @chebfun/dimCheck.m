function out = dimCheck(f, g, ~)
%DIMCHECK   Check dimension compatability of two CHEBFUN objects.
%   In MATLAB 2015b and below DIMCHECK(F, G) returns:
%       1 if numColumns(F) == numColumns(G)
%       0 otherwise.
%
%   In MATLAB 2016a and above DIMCHECK(F, G) returns:
%       1 if numColumns(F) == numColumns(G)
%       2 if numColumns(F) == 1
%       3 if numColumns(G) == 1
%       0 otherwise.
%
%   DIMCHECK(F, G, 1) will thrown a dimension error rather than return 0. 
%
% See also NUMCOLUMNS, SIZE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

nf = numColumns(f);
ng = numColumns(g);

if ( nf == ng )
    out = 1;
    
elseif ( nf == 1 )
    out = 2;
    
elseif ( ng == 1 )
    out = 3;
    
else
    out = 0;
    
end

if ( (out > 0) && verLessThan('matlab', '8.6') )
   out = 0; 
end

if ( nargin == 3 && out == 0 )
    error('CHEBFUN:CHEBFUN:dimCheck:dims', 'Matrix dimensions must agree.');
end

end