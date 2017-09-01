function varargout = rank(F)
%RANK   Trilinear or Tucker rank of a CHEBFUN3 object.
%   [rX, rY, rZ] = RANK(F) returns the size of the core tensor, i.e., rX is the
%   number of columns in factor quasimatrix F.cols, rY is the number of 
%   columns in F.rows and rZ is the number of columns in F.tubes.
%
%   r = RANK(F) returns the maximum entry in the size of the core tensor.
%
%   An exception is the zero CHEBFUN3 object for which rank is defined to 
%   be zero. 
%
% See also CHEBFUN3/HOSVD.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

r = size(F.core);
if ( numel(r) < 3 )
    % Developer note: If the input function handle is bivariate, then the 
    % core tensor reduces to be a 2D matrix. Try e.g. 
    % >> f = chebfun3(@(x, y, z) cos(x + y))
    % If this has happened, then manually force it to be a tensor.
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