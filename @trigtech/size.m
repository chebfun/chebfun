function [size1, size2] = size(f, varargin)
%SIZE   Size of a TRIGTECH.
%   [S1, S2] = SIZE(F) returns S1, the number of values at equally spaced
%   points used to define F, and S2, the number of columns in F.
%
%   S = SIZE(F) returns the same as above in a 1x2 vector, S = [S1, S2].
%
% See also LENGTH.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% The size of a TRIGTECH is the size of its vector of values.
size1 = size(f.values, varargin{:});

% Return two outputs:
if ( nargout == 2 )
    size2 = size1(2);
    size1 = size1(1);
end

end