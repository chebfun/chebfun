function varargout = solvebvp(varargin)
%SOLVEBVP   Solve CHEBOP BVP system.
%   SOLVEBVP is a wrapper for CHEBOP/MLDIVIDE and is maintained for backwards
%   compatability. See CHEBOP/MLDIVIDE documentation for further details.
%
% See also CHEBOP/MLDIVIDE.

% Call CHEBOP/MLDIVIDE:
[varargout{1:nargout}] = mldivide(varargin{:});

end