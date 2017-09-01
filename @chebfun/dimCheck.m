function out = dimCheck(f, g, ~)
%DIMCHECK   Check dimension compatability of two CHEBFUN objects.
%   In MATLAB 2016a and below DIMCHECK(F, G) returns:
%       1    if numColumns(F) == numColumns(G)
%     error  otherwise.
%
%   In MATLAB 2016b and above DIMCHECK(F, G) returns:
%       1    if numColumns(F) == numColumns(G)
%      -1    if numColumns(F) == 1 or numColumns(G) == 1
%     error  otherwise.
%
%   DIMCHECK(F, G, 1) will return 0 rather than throw a dimension error.
%
% See also NUMCOLUMNS, SIZE.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

nf = numColumns(f);
ng = numColumns(g);

% Check dimensions:
if ( nf == ng )
    out = 1;
elseif ( nf == 1 || ( ng == 1 ) )
    out = -1;
else
    out = 0;
end

% Adjust for MATLAB version if out = -1:
if ( (out == -1 ) && verLessThan('matlab', '9.1') )
   out = 0; 
end

if ( out == 0 && nargin == 2 )
    % Throw an error if requested (as the caller):
    ME = MException('CHEBFUN:CHEBFUN:dimCheck:dim', ...
        'Matrix dimensions must agree.');
    throwAsCaller(ME)
end

end