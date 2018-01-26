function varargout = quiver(N, varargin)
% QUIVER   Draw phase plot diagrams, based on ODEs specified with CHEBOPS
%
% Calling sequence:
%   H = QUIVER(N, AXIS, 'OPT1', VAL1, ...)
%
% Here, the inputs are:
%   N    : A chebop, whose N.op arguments specifies a first or second order
%          scalar ODE, or a coupled system of two first order ODEs.
%   AXIS : A 4-vector with elements [XMIN XMAX YMIN YMAX] that specify the
%          rectangular region shown on the phase plot. If none is passed, the
%          default values [-1 1 -1 1] are used.
%
% It is possible to pass the method various option pairs of the form
% 'OPTIONNAME', OPTIONVALUE. The options supported are:
%   'XPTS'      : An integer, specifying the resolution of the x-axis for the
%                 quiver plot. Default value: 10.
%   'YPTS'      : An integer, specifying the resolution of the y-axis for the
%                 quiver plot. Default value: 10.
%   'NORMALIZE' : A Boolean, which determines whether the arrows on the quiver
%                 plot or normalized all to have the same length. Default value:
%                 false.
%   'SCALE':      By default, quiver automatically scales the arrows to fit
%                 within the grid. By passing a SCALE argument S, the arrows are
%                 fitted within the grid and then stretched by S. Use S = 0 to
%                 plot the arrows without the automatic scaling. See the
%                 documentation for the built-in MATLAB QUIVER for more
%                 information.
%   LINESPEC      Specifies options for the quiver plot. The LINESPEC argument
%                 can either be a string supported by the MATLAB PLOT command,
%                 such as 'ro', or a parameter/value pair to specify additional
%                 properties of the quiver lines, such as 
%                   quiver(N,[0 2 0 2], 'k', 'linewidth', 1.4).
%
% The optional output is
%   H   : A quivergroup handle.
%
% Note: The CHEBOP QUIVER command works by reformulating higher order problems
% as coupled first order systems, evaluating the resulting first order system at
% grid that should be interpreted as values of u and u', then calling the
% built-ing MATLAB QUIVER method on the results. In the case of first order
% scalar problems, the grid should be interpreted as values of t and u.
%
% Example 1 -- van der Pol equation (second order ODE)
%   N = chebop(0,100);
%   N.op = @(t,u) diff(u, 2) - 3*(1-u^2)*diff(u) + u;
%   quiver(N,[-2.5 2.5 -5.5 5.5], 'r', 'xpts', 40, 'ypts', 40, 'scale', .5, ...
%       'normalize', true, 'linewidth', 1.5)
%   hold on % Plot a particular solution on top of quiver plot
%   N.lbc = [2; 0];
%   u = N\0;
%   plot(u, diff(u),'linewidth',2)
%
% Example 2 -- Lotka-Volterra (first order coupled system)
%   N = chebop(@(t,u,v) [diff(u)-2.*u+u.*v; diff(v)+v-u.*v], [0 4]);
%   quiver(N, [0 2.5 0 4])
%   hold on
%   N.lbc = @(u,v) [u - 0.5; v - 1]; % Initial populations
%   [u, v] = N\0;
%   plot(u, v, 'linewidth', 2)
%   plot(0.5, 1,'m*','markersize',15) % Mark initial condition
%
% Example 3 -- Slopefield for a first order problem
%   N = chebop(@(t,u) diff(u)-sin(t)*u);
%   quiver(N,[-1.2*pi 1.2*pi -1 1])


% Copyright 2017 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

% Set default values:
scale = 1;
normalize = false;
u0 = [];
xpts = 20;
ypts = 20;
axisLims = [-1 1 -1 1];

% The limits for the axes have to appear in the first entry of varargin if
% they're specified. They'll always be a 4-vector.
if ( ~isempty(varargin) && length(varargin{1}) == 4 )
    axisLims = varargin{1};
    % Throw away the axisLims:
    varargin(1) = [];
end

% Cell for linespecs
linespec = {};

% Parse VARARGIN, go through all elements:
while ( ~isempty(varargin) )    
    if ( ~ischar(varargin{1}) && ~isnumeric(varargin{2}) )
        error('followpath:inputArgument','Incorrect options input arguments');
    end
    
    throwAwayLength = 2;
    switch lower(varargin{1})
        case 'normalize'
            normalize = varargin{2};
        case 'xpts'
            xpts = varargin{2};
        case 'ypts'
            ypts = varargin{2};
        case 'scale'
            scale = varargin{2};
        otherwise
            % Must have gotten something to do with linespec
            findNorm = strcmpi(varargin, 'normalize');
            findXpts = strcmpi(varargin, 'xpts');
            findYpts = strcmpi(varargin, 'ypts');
            findScale = strcmpi(varargin, 'scale');
            optionPos = findNorm | findXpts | findYpts | findScale;
            nextOption = find(optionPos, 1, 'first');
            % If nextOption is empty, we only got linespec options left
            if ( isempty(nextOption) )
                nextOption = length(varargin) + 1;
            end
            linespec = [linespec, varargin{1:nextOption - 1}];
            throwAwayLength = nextOption - 1;
    end
    
    % Throw away option name and argument and move on:
    varargin(1:throwAwayLength) = [];
end

% Extract the x and y limits:
xl = axisLims(1:2);
yl = axisLims(3:4);

% If ylim is empty, we solve the problem to obtain a range for plotting on:
if ( isempty(xl) )
    u0 = N\0;
    xl = 1.1*minandmax(u0);
end

if ( isempty(yl) )
    if ( isempty(u0) )
        u0 = N\0;
    end
    yl = 1.1*minandmax(diff(u0));
end

% Convert the operator in N to first order.
[firstOrderFun, ~, ~, ~, diffOrd] = treeVar.toFirstOrder(N.op, 0, N.domain);

% Vectors for constructing a meshgrid:
y1 = linspace(xl(1), xl(end), xpts);
y2 = linspace(yl(1), yl(end), ypts);

% Get a meshgrid for points of interests in the phase plane.
[x,y] = meshgrid(y1, y2);
u = zeros(size(x));
v = zeros(size(x));

% Phase plane portraits really only make sense for autonomous systems, which
% shouldn't depend on t, hence, we simply take t = 0 for evaluating.
t = 0;

% Check if we got passed a system of too high order:
numEquations = length(strfind(func2str(firstOrderFun),';'));

if ( numEquations > 2 )
    error('CHEBFUN:CHEBOP:quiver:tooHighOrder', ...
        ['The ODE passed to chebop/quiver must either be a scalar, second ' ...
        'order ODE or a system of two first order equations.'])
end

% Are we plotting a slope field (cf. #2238), or a standard phase plane?
if( (length(diffOrd) == 1) && (diffOrd == 1))   % Slope field
    for i = 1:numel(x)
        res = firstOrderFun(x(i), y(i));
        u(i) = 1;
        v(i) = res(1);
    end
else                                            % Phase plane
    for i = 1:numel(x)
        res = firstOrderFun(t,[x(i); y(i)]);
        u(i) = res(1);
        v(i) = res(2);
    end
end

if ( normalize )
    % Make all arrows equal length:
    nrm = sqrt(u.^2 + v.^2);
    u = u./nrm;
    v = v./nrm;
end

h = quiver(x, y, u, v, scale,linespec{:});

% Set x and y limits:
xlim(xl);
ylim(yl);
if ( nargout > 0 )
    varargout{1} = h;
else
    
end
