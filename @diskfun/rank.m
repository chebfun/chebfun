function varargout = rank(varargin)
%RANK   Rank of a DISKFUN.
%   RANK(F) produces an estimate of the rank of F.
%
%   RANK(F, TOL) is the number of singular values of F greater than TOL*N,
%   where N is the first singular value of F.
%
% See also DISKFUN/LENGTH.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

[varargout{1:nargout}] = rank@separableApprox(varargin{:});

end