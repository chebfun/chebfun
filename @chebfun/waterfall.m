function varargout = waterfall(varargin)
%WATERFALL   Waterfall plot for CHEBFUN object.
%   WATERFALL(U), or WATERFALL(U, T) where LENGTH(T) = MIN(SIZE(U)), plots a
%   "waterfall" plot of an array-valued CHEBFUN or quasimatrix. Unlike the
%   standard Matlab WATERFALL method, CHEBFUN/WATERFALL does not fill in the
%   column planes with opaque whitespace or connect edges to zero. Instead,
%   horizontal slices are connected by a semi-transparent egde.
%
%   Additional plotting options can also be passed, for example WATERFALL(U, T,
%   'linewidth', 2). Additional options include 'EdgeColor', 'EdgeAlpha',
%   'FaceAlpha', and 'FaceColor'. See the built-in WATERFALL method for details.
%
% See also PLOT, PLOT3.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% % First input might be a figure handle:
% [cax, varargin] = axescheck(varargin{:});
% if ( ~isempty(cax) )
%     axes(cax);
% end

% First input is now the CHEBFUN:
u = varargin{1};
varargin(1) = [];

% Some defaults:
faceAlpha = .75;
faceColor = 'w';

% Parse inputs:
k = 1;
while ( k <= numel(varargin) )
    if ( strcmpi(varargin{k}, 'numpts') )
        % This no longer does anything.
%         numpts = varargin{k+1};
        varargin(k:k+1) = [];
    elseif ( strcmpi(varargin{k}, 'simple') )
        % Note, this no longer does any thing. Simple is always true.
        varargin(k) = [];
    elseif ( strcmpi(varargin{k}, 'facealpha') )
        faceAlpha = varargin{k+1};
        varargin(k:k+1) = [];
    elseif ( strcmpi(varargin{k}, 'facecolor') )
        faceColor = varargin{k+1};        
        varargin(k:k+1) = [];
    elseif ( strcmpi(varargin{k}, 'fill') )
        varargin(k) = [];
        warning('CHEBFUN:CHEBFUN:waterfall:fill', ...
            'The ''fill'' option is no longer supported.');
        warning('off', 'CHEBFUN:CHEBFUN:waterfall:fill');
    else
        k = k + 1;
    end
end
nargin = length(varargin) + 1;

% Convert quasimatrix to an array-valued CHEBFUN.
u = quasi2cheb(u);    

% Deal with column vectors:
isTrans = u(1).isTransposed;
if ( isTrans )
    u = u.';
end
n = min(size(u, 2));

% Grab t if it is given, if not, default to linear:
if ( (nargin > 1) && (isnumeric(varargin{1})) && ...
        (length(varargin{1}) == size(u, 2)) )
    t = varargin{1}; 
    t = t(:).';
    varargin(1) = [];
else
    t = 1:n;
end

% Some error checking:
if ( length(t) ~= n )
    error('CHEBFUN:CHEBFUN:waterfall:sizeT', ...
        'Length of T should equal the number of quasimatrices in U');
end
if ( ~isreal(u) || ~all(isreal(t)) )
    warning('CHEBFUN:CHEBFUN:waterfall:imaginary',...
        'Imaginary parts of complex T and/or U arguments ignored');
    u = real(u); 
    t = real(t);
end    

% Get the data:
data = plotData(u);
uu = data.yLine;
% Repeat x and t for each column:
xx = repmat(data.xLine, 1, n);
tt = repmat(t, length(xx(:,1)), 1);
% Mask the NaNs:
mm = find(isnan(uu));
uu(mm) = uu(mm + 1);
uu(1) = uu(1) + eps;

% If we're given a scalar, deal with it: 
if ( min(size(uu)) == 1 )
    xx = [xx, xx]; 
    tt = [tt, tt];
    uu = [uu, uu]; 
    t =  [t, t+1];
end

% Store the state of HOLD:
holdState = ishold;

% Plot the whitespace:
h1 = mesh(xx.', tt.', uu.', ...
    'EdgeColor', 'none', 'Facealpha', faceAlpha, 'FaceColor', faceColor); 
hold on

% Intersperse with data NaNs to prevent fill:
xx = repmat(xx, 1, 4);
mid = repmat((t(1:end-1) + t(2:end))/2, size(tt, 1), 1);
tt = reshape([tt                              ; 
              tt+eps                          ; 
              [mid NaN*tt(:,end)]             ; 
              [tt(:,2:end)-eps NaN*tt(:,end)]], size(xx));
uu = reshape([uu                      ;    
              uu                      ; 
              NaN*uu                  ; 
              [uu(:,2:end) uu(:,end)]], size(xx));

% Plot the edges:
h2 = mesh(xx.', tt.', uu.', varargin{:});

% Reset the HOLD state:
if ( ~holdState )
    hold off
end

% Output figure handles:
if ( nargout > 0 )
    h = [h1 ; h2];
    varargout{1} = h;
end

end
