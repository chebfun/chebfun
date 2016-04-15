function f = cumsum3(f)
%CUMSUM2   Triple indefinite integral of a CHEBFUN3.
%   F = CUMSUM3(F) returns the triple indefinite integral of a CHEBFUN3. 
% That is
%                 z  y  x
%                /  /  /
%   CUMSUM3(F) =|  |  |   F(x,y,z) dx dy dz
%              /  /  /
%             e  c  a
%
%   for  (x,y,z) in [a,b] x [c,d] x [e,g], where 
%   [a,b] x [c,d] x [e,g] is the domain of F.
% 
%   See also chebfun3/cumsum, chebfun3/cumsum2, chebfun3/sum, 
%   chebfun3/sum and chebfun3/sum3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Check for empty:
if ( isempty(f) ) 
    f = [];
    return
end

f.cols = cumsum(f.cols);     % CUMSUM along the cols.
f.rows = cumsum(f.rows);     % CUMSUM along the rows.
f.tubes = cumsum(f.tubes);   % CUMSUM along the tubes.

end