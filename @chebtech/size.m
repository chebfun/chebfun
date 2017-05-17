function [size1, size2] = size(f, varargin)
%SIZE   Size of a CHEBTECH.
%   [S1, S2] = SIZE(F) returns S1, the number of values at Chebyshev points used
%   to define F, and S2, the number of columns in F.
%
%   S = SIZE(F) returns the same as above in a 1x2 vector, S = [S1, S2].
%
% See also LENGTH.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% The size of a CHEBTECH is the size of its vector of coeffs.
size1 = size(f.coeffs, varargin{:});

% Return two outputs:
if ( nargout == 2 )
    size2 = size1(2);
    size1 = size1(1);
end

end
