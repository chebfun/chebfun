function h = quiver(N, varargin)

% Set default values
rhs = 0;
xl = [];
yl = [];
normalize = false;
u0 = [];
xpts = 10;
ypts = 10;


% Parse VARARGIN
while ~isempty(varargin)    % Go through all elements
    if ~ischar(varargin{1}) && ~isnumeric(varargin{2})
        error('followpath:inputArgument','Incorrect options input arguments');
    end
    val = varargin{2};
    switch lower(varargin{1})
        case 'rhs'
            rhs = val;
        case 'xlim'
            xl = val;
        case 'ylim'
            yl = val;
        case 'normalize'
            normalize = val;
        case 'xpts'
            xpts = val;
        case 'ypts'
            ypts = val;
    end
    
    % Throw away option name and argument and move on
    varargin(1:2) = [];
end

% If ylim is empty, we solve the problem to obtain a range for plotting on
if ( isempty(xl) )
    u0 = N\rhs;
    xl = 1.1*minandmax(u0);
end

if ( isempty(yl) )
    if ( isempty(u0) )
        u0 = N\rhs;
    end
    yl = 1.1*minandmax(diff(u0));
end

% Convert the operator in N to first order
firstOrderFun = treeVar.toFirstOrder(N.op, rhs, N.domain);

%%

y1 = linspace(xl(1), xl(end), xpts);
y2 = linspace(yl(1), yl(end), ypts);

% Get a meshgrid for points of interests in the phase plane.
[x,y] = meshgrid(y1, y2);
u = zeros(size(x));
v = zeros(size(x));

% Phase plane portraits really only make sense for autonomous systems, which
% shouldn't depend on t, hence, we simply take t = 0 for evaluating
t=0; 

% TODO: Could probably vectorize this, with reshapes. For now, just loop.
for i = 1:numel(x)
    res = firstOrderFun(t,[x(i); y(i)]);
    u(i) = res(1);
    v(i) = res(2);
end

if (normalize)
    % Make all arrows equal length
    nrm = sqrt(u.^2 + v.^2);
    u = u./nrm;
    v = v./nrm;
end

h = quiver(x, y, u, v);
xlim(xl);
ylim(yl);
end
