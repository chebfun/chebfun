function varargout = transpose(varargin) %#ok<STOUT,*INUSD>
%TRANSPOSE   DELTAFUN objects are not transposable.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

error('DELTAFUN:transpose:notpossible', ...
    'DELTAFUN objects are not transposable.')

end
