function varargout = rank(varargin)
%RANK      Rank of a CHEBFUN2.
%   RANK(F) produces an estimate of the rank of the approximant F.
%
%   RANK(F, TOL) is the number of singular values of F greater than TOL/N, where
%   N is the first singular value of F.
%
% See also LENGTH.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = rank@separableApprox(varargin{:});

end
