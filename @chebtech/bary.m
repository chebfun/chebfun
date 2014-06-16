function fx = bary(x, fvals, xk, vk)
%BARY   Barycentric interpolation formula.
%   BARY(X, FVALS, XK, VK) uses the 2nd form barycentric formula with weights VK
%   to evaluate an interpolant of the data {XK, FVALS(:,k)} at the points X.
%   Note that X, XK, and VK should be column vectors, and FVALS, XK, and VK
%   should have the same length.
%
% See also CHEBTECH.CHEBPTS, CHEBTECH.BARYWTS, CHEBTECH.FEVAL, CHEBTECH.CLENSHAW.

% [TODO]: Move this to either the trunk/ folder or an @utils class.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs
[n, m] = size(fvals);
sizex = size(x);
ndimsx = ndims(x);

if ( ~all(sizex) )
    fx = x;
    return
end

% Check that input is a column vector:
if ( (ndimsx > 2) || (sizex(2) > 1) ) 
    warning('CHEBFUN:CHEBTECH:bary:colvec', 'Input should be a column vector.');
    x = x(:);
end

% The function is a constant.
if ( n == 1 )
    fx = repmat(fvals, length(x), 1);
    return
end

% The function is NaN.
if ( any(isnan(fvals)) )
    fx = NaN(length(x), m);
    return
end

% The main loop:
if ( numel(x) < 4*length(xk) )  % Loop over evaluation points
    % Note: The value "4" here was detemined experimentally.
    
    % Initialise return value:
    fx = zeros(size(x, 1), m);

    % Loop:
    for j = 1:numel(x),
        xx = vk ./ (x(j) - xk);
        fx(j,:) = (xx.'*fvals) / sum(xx);
    end
else                            % Loop over barycentric nodes
    % Initialise:
    num = zeros(size(x, 1), m);
    denom = num;

    % Loop:
    for j = 1:length(xk),
        tmp = (vk(j) ./ (x - xk(j)));
        num = num + bsxfun(@times, tmp, fvals(j,:));
        denom = bsxfun(@plus, denom, tmp);
    end
    fx = num ./ denom;
end

% Try to clean up NaNs:
for k = find(isnan(fx(:,1)))'       % (Transpose as Matlab loops over columns)
    index = find(x(k) == xk, 1);    % Find the corresponding node
    if ( ~isempty(index) )
        fx(k,:) = fvals(index,:);   % Set entry/entries to the correct value
    end
end

end
