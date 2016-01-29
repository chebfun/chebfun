function [y, x] = min(f, flag, dim)
%MIN   Minimum values of a CHEBFUN.
%   MIN(F) and MIN(F, 'global') return the minimum value of the CHEBFUN F.
%
%   [Y, X] = MIN(F) returns also a point X such that F(X) = Y.
%
%   [Y, X] = MIN(F, 'local') returns not just the global minimum value but all
%   of the local minima.
%
%   If F is complex-valued, absolute values are taken to determine the minima,
%   but the resulting values correspond to those of the original function.
%
%   If F is array-valued, then the columns of X and Y correspond to the columns
%   of F. NaNs are used to pad Y and X when the 'local' flag is used and the
%   columns are not of the same length.
%
%   H = MIN(F, G), where F and G are CHEBFUNs defined on the same domain,
%   returns a CHEBFUN H such that H(x) = min(F(x), G(x)) for all x in the
%   domain of F and G. Alternatively, either F or G may be a scalar.
%
%   MAX(F, [], DIM) computes the maximum of the CHEBFUN F in the dimension DIM.
%   If DIM = 1 and F is a column CHEBFUN or DIM = 2 and F is a row CHEBFUN, this
%   is equivalent to MAX(F). Otherwise, MAX(F, [], DIM) returns a CHEBFUN which
%   is the maximum across the discrete dimension of F. For example, if F is a
%   quasimatrix with two columns, MAX(F, [], 2) = MAX(F(:,1), F(:,2)).
%
% See also MAX, MINANDMAX, ROOTS.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Trivial empty case.
if ( isempty(f) )
    x = [];
    y = [];
    return
end

if ( nargin == 3 )
    if ( ~any(dim == [1 2]) )
        error('CHEBFUN:CHEBFUN:min:badDim', ...
            'DIM input to CHEBFUN MIN must be 1 or 2.');
    end

    if ( dim ~= 1 + f(1).isTransposed )
        % Take min across discrete dimension of a quasimatrix:
        f = cheb2cell(f);
        y = realmax;
        for k = 1:numel(f)
            y = minOfTwoChebfuns(y, f{k});
        end
        y = merge(y);
        return
    end
end

if ( (nargin > 2) && isempty(flag) )
    % MAX(F, [], 1).
    flag = 'global';
end

if ( (nargin == 1) || strcmp(flag, 'global') ) 
    % MIN(F) or MIN(F, 'global')
    [y, x] = globalMin(f);    
    
elseif ( isa(flag, 'chebfun') || isnumeric(flag) )
    % MIN(F, G)
    y = minOfTwoChebfuns(f, flag);
    
elseif ( strcmp(flag, 'local') )
    % MIN(F, 'local')
    [y, x] = localMin(f);
    
else
    error('CHEBFUN:CHEBFUN:min:flag', 'Unrecognized flag.');
    
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% GLOBAL MINIMUM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y, x] = globalMin(f)

% Call MINANDMAX():
[y, x] = minandmax(f);

% Extract the minimum:
y = y(1,:);
x = x(1,:);

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% LOCAL MINIMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y, x] = localMin(f)

% Call MINANDMAX():
[y, x] = minandmax(f, 'local');

% Determine which are minima.

ends = f(1).domain([1, end]).'; % Endpoints of the domain are special.
f = mat2cell(f); % Convert f into a cell of scalar-valued CHEBFUNs.

% Loop over the columns:
for k = 1:numel(f)
    % Compute 1st and 2nd derivatives.
    dfk = diff(f{k});
    dfk2 = diff(dfk);

    % For interior extrema, look at 2nd derivative:
    minimaLoc = feval(diff(f{k}, 2), x(:,k)) > 0;

    % For end-points, look at 1st derivative:
    dfk_ends = feval(dfk, ends);
    endptMinLoc = dfk_ends.*[1, -1]' > 0;

    % If 1st derivative is small at an endpoint, assume it's zero and try to
    % use 2nd derivative to determine if it's a minimum.
    %
    % [TODO]:  What if the 2nd derivative is zero, so that rounding error
    % precludes us from accurately determining the sign?
    smallEndDer = abs(dfk_ends) < 1e3*vscale(dfk)*eps;
    endptMinLoc(smallEndDer) = feval(dfk2, ends(smallEndDer)) > 0;

    minimaLoc(1) = endptMinLoc(1);
    minimaLoc(x(:,k) == ends(2)) = endptMinLoc(2);

    % Set points corresponding to local maxima to NaN:
    y(~minimaLoc,k) = NaN;
    x(~minimaLoc,k) = NaN;

    % Sort the result
    [x(:,k), minimaLoc] = sort(x(:,k));
    y(:,k) = y(minimaLoc,k);
end

% Remove any rows which contain only NaNs.
x(all(isnan(x), 2),:) = []; 
y(all(isnan(y), 2),:) = [];

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% MIN(F, G) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = minOfTwoChebfuns(f, g)
% Return the function h(x) = max(f(x), g(x)) for all x. 

% If one is complex, use abs(f) and abs(g) to determine which function values to
% keep. (experimental feature)
if ( isreal(f) && isreal(g) && (nargin < 3) )
	S = sign(f - g);
else
	S = sign(abs(f) - abs(g));
end

% Heaviside function (0 where f > g, 1 where f < g);
H = 0.5*(S + 1);
notH = 0.5*(1 - S); % ~H.

% Combine for output:
h = notH.*f + H.*g;

% [TODO]: Enforce continuity?

% TODO: Why simplify?
% Simplify:
h = simplify(h);

end
