function varargout = sampleTest(varargin)
%SAMPLETEST   Test an evaluation of input OP against a SPHEREFUN.
%
%   SAMPLETEST(F, SAMPLEOP, TOL) evaluates both the function OP and it
%   SPHEREFUN representation F at several points in it's domain. The difference of
%   these values is computed, and if this is sufficiently small the test 
%   passes and returns TRUE. If the difference is large, it returns FALSE.
% 
%   SAMPLETEST(F, SAMPLEOP, TOL, FLAG) is the same as above if FLAG = 0. 
%   If FLAG = 1 then the OP is assumed to be unvectorized. 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = sampleTest@separableApprox(varargin{:});

end
