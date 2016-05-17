function varargout = sampleTest(varargin)
%SAMPLETEST   Test an evaluation of input OP against a CHEBFUN2.
%
%   SAMPLETEST(F, SAMPLEOP, TOL) evaluates both the function OP and its
%   CHEBFUN2 representation F at several points in its domain. The difference of
%   these values is computed, and if this is sufficiently small the test
%   passes and returns TRUE. If the difference is large, it returns FALSE.
%
%   SAMPLETEST(F, SAMPLEOP, TOL, FLAG) is the same as above if FLAG = 0.
%   If FLAG = 1 then the OP is assumed to be unvectorized.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = sampleTest@separableApprox(varargin{:});

end
