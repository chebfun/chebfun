function f = permute(f, order)
%PERMUTE   Permute CHEBFUN array dimensions.
%   G = permute(F, ORDER) rearranges the dimensions of A so that they are in the
%   order specified by the vector ORDER. The array produced has the same values
%   as A but the order of the subscripts needed to access any particular element
%   are rearranged as specified by ORDER. Since CHEBFUN objects ony have two
%   dimensions, ORDER must be one of [1, 2, ...] or [2, 1, ...]. In the first
%   case, G = F, and in the second, G = F.';
%
% See also TRANSPOSE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% For simplicity, we attempt to permute a small matrix and catch any errors
% in inputs this way.
tmp = eye(2);
permute(tmp, order);

% If we got here, then v = [1, 2] or [2, 1].
if ( order(1) == 2 )
    % Transpose f:
    f = f.';
end

end
