function [size1, size2] = size(f,varargin)
%LENGTH	  Length of a FUN.
%   SIZE(F), where F is a FUN, is the size of the ONEFUN of F.
%
%   See also LENGTH.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% The size of a FUN is the size of its onefun.

if nargout == 1
    size1= size(f.onefun, varargin{:});
else
    [size1, size2] = size(f.onefun, varargin{:});
end

end
