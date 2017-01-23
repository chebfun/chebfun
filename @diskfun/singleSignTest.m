function varargout = singleSignTest(varargin)
%SINGLESIGNTEST   Heuristic check to see if a DISKFUN changes sign. 
%   SINGLESIGNTEST(F) returns 1 if the values of F on a tensor grid are of 
%   the same sign.
%
%   [OUT, WZERO] = SINGLESIGNTEST( F ), returns WZERO = 1 if a zero has 
%   been found.
%
%   The algorithm works by sampling F on a tensor-grid and checking if
%   those values are of the same sign. This command is mainly for internal 
%   use in DISKFUN commands.
%
% See also DISKFUN/ABS, DISKFUN/SQRT, and DISKFUN/LOG.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = singleSignTest@separableApprox(varargin{:});

end