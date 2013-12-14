function [size1, size2] = size(f, varargin)
%SIZE   Size of a DELTAFUN.
%   [S1, S2] = SIZE(F) returns the size of the funPart of F.
%
%   S = SIZE(F) returns the same as above in a 1x2 vector, S = [S1, S2].
%
% See also LENGTH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% The size of a DELTAFUN is the size of its funPart.
size1 = size(f.funPart, varargin{:});

% Return two outputs:
if ( nargout == 2 )
    size2 = size1(2);
    size1 = size1(1);
end

end
