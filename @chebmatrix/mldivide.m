function varargout = mldivide(A, b, varargin)

% TODO: Document this.

[varargout{1:nargout}] = mldivide(linop(A), b, varargin{:});

end