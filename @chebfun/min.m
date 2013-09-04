function [y, x] = min(f, flag)
%MIN    Minimum values of a CHEBFUN.
%   MIN(F) returns the minimum value of the chebfun F.
%
%   [Y, X] = MIN(F) returns also point X such that F(X) = Y.
%
%   [Y, X] = MIN(F, 'local') returns not just the global minimum value, but all
%   of the local minima.
%
%   If F is complex-valued, absolute values are taken to determine the minima,
%   but the resulting values correspond to those of the original function.
%
%   If F is array-valued, then the columns of X and Y correspond to the columns
%   of F. NaNs are used to pad Y and X when the 'local' flag is used and the
%   columns are not of the same length
%
% See also MAX, MINANDMAX, ROOTS.
%
% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.


if ( nargin == 1 )
    [y, x] = globalMin(f);

elseif ( isa(flag, 'chebfun') )
    % [TODO]: Implement this. (Requires SIGN())
    y = minOfTwoChebfuns(f, flag);
 
else
    [y, x] = localMin(f);

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
[y, x] = minandmax(f, flag);

% Determine which are minima.

ends = f.domain([1, end]).'; % Endpoints of the domain are special.
f = mat2cell(f); % Convert f into a cell of scalar-valued CHEBFUNs.

% Loop over the FUNs:
for k = 1:numel(f)
    % For interior extrema, look at 2nd derivative:
    idx = feval(diff(f{k}, 2), x(:,k)) > 0;
    % For end-points, look at 1st derivative:
    idx2 = feval(diff(f{k}), ends).*[1, -1]' > 0;
    idx(1) = idx2(1);
    idx(x(:,k) == ends(2)) = idx2(2);
    % Set points corresponding to local maxima to NaN:
    y(~idx,k) = NaN;
    x(~idx,k) = NaN;
    % Sort the result
    [x(:,k), idx] = sort(x(:,k));
    y(:,k) = y(idx,k);
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
if ( isreal(f) && isreal(g) && nargin < 3 )
	S = sign(f - g);
else
	S = sign(abs(f) - abs(g));
end

% Heaviside function (0 where f > g, 1 where f < g);
H = ((S+1)/2);
notH = ((1-S)/2); % ~H.

% Combine for output:
h = notH.*f + H.*g;

% [TODO]: Enforce continuity?

end