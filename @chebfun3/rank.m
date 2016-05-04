function varargout = rank(F)
%RANK   Multilinear rank of a CHEBFUN3 object F.
%
%   If three outputs are asked for, then it is the size of the core tensor, 
%   i.e., r = [size(F.cols,2), size(F.rows,2), size(F.tubes,2)];
%   If just one output is asked for, then r is the max size of the core 
%   tensor. 
%
%   See also CHEBFUN3/HOSVD.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


r = size(F.core);
if numel(r)<3
    r = [r, 1];
end

% Check for zero function.
if ( all(F.core(:) == 0) )
    r = [0, 0, 0]; 
end

% Output:
if ( nargout <= 1 )
    varargout = {max(r(:))};
else
    varargout = {r(1), r(2), r(3)};
end

end
