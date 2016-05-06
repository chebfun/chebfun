function varargout = surfc(u, varargin)
%SURFC   Combination surf/contour plot for CHEBFUN objects.
%   SURFC(...) is the same as SURF(...) except that a contour plot is drawn
%   beneath the surface.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

transState = u(1).isTransposed;
if ( transState )
    u = u.';
end
n = numColumns(u);
t = 1:n;

if ( (nargin > 1) && isnumeric(varargin{1}) && (length(varargin{1}) == size(u,2)) )
    t = varargin{1};
    t = t(:).';
    varargin(1) = [];
end

if ( (numel(varargin) > 1) && strcmpi(varargin{1}, 'numpts') )
    numpts = varargin{2};
    varargin(1:2) = [];
    warning('CHEBFUN:CHEBFUN:surf:numpts', 'NUMPTS option is deprecated.');
end

if ( length(t) ~= n )
    error('CHEBFUN:CHEBFUN:surfc:szet', ...
        'Length of T should equal the number of quasimatrices in U');
end

if ( ~isreal(u) || ~all(isreal(t)) )
    warning('CHEBFUN:CHEBFUN:surfc:imaginary', ...
        'Imaginary parts of complex T and/or U arguments ignored');
    u = real(u); t = real(t);
end

if ( n == 1 )
    [varargout{1:nargout}] = waterfall(u, t, varargin{:});
    return
end

% Convert to a quasimatrix:
u = quasi2cheb(u);
if ( numel(u) > 1 )
    error('CHEBFUN:CHEBFUN:surfc:quasi', ...
        'SURF does not support quasimatrices.');
end

% Get the data:
data = plotData(u);
uu = data.yLine;
xx = repmat(data.xLine, 1, n);
tt = repmat(t, length(xx(:,1)), 1);

defaults = {'edgealpha', 0};

% Plot the surface:
if ( ~transState )
    h = surfc(xx.', tt.', uu.', defaults{:}, varargin{:});
else
    h = surfc(xx.', tt.', uu, defaults{:}, varargin{:});
end
shading interp

if ( nargout > 0 )
    varargout{1} = h;
end

end
