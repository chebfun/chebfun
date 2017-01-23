function varargout = size(f, varargin)
%SIZE   Size of a DELTAFUN.
%   [S1, S2] = SIZE(F) returns the size of the funPart of F.
%
%   S = SIZE(F) returns the same as above in a 1x2 vector, S = [S1, S2].
%
% See also LENGTH.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% The size of a DELTAFUN is the size of its funPart.
varargout{1:nargout} = size(f.funPart, varargin{:});

end
