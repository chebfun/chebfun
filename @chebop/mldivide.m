function varargout = mldivide(varargin)
%\    MLDIVIDE   Solve CHEBOP BVP system.
%   MLDIVIDE is a convenient wrapper for CHEBOP/SOLVEBVP, but is limited in that
%   it only supports a single output. See CHEBOP/SOLVEBVP documentation for
%   further details.
%
% See also CHEBOP/SOLVEBVP.

% Call CHEBOP/MLDIVIDE:
[varargout{1:nargout}] = solvebvp(varargin{:});

end