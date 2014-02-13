function varargout = domain(f, varargin)
% DOMAIN(F) where F is an ADCHEBFUN is the same as DOMAIN(F.FUNC)
[varargout{1:nargout}] = domain(f.func, varargin{:});
end
