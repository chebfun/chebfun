function [size1, size2] = size(f, varargin)
%SIZE   Size of a SINGFUN.
%   [S1, S2] = SIZE(F) returns the size of the smoothPart of F.
%
%   S = SIZE(F) returns the same as above in a 1x2 vector, S = [S1, S2].
%
% See also LENGTH.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% The size of a SINGFUN is the size of its smooth part:
size1 = size(f.smoothPart, varargin{:});

% Return two outputs:
if ( nargout == 2 )
    size2 = size1(2);
    size1 = size1(1);
end

end
