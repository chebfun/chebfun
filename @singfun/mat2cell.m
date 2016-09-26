function g = mat2cell(f, M, N)
%MAT2CELL   SINGFUN MAT2CELL.
%   As SINGFUN does not support array-valued objects, G = MAT2CELL(F, M, N) 
%   returns G = {F}, regardless of M and N, to make it compatible with the
%   prototype for MAT2CELL established by ONEFUN.
%
% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = {f};
    
end
