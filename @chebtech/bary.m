function fx = bary(x, fvals, xk, vk)
%BARY   Barycentric interpolation formula.
%   BARY(X, FVALS, XK, VK) uses the 2nd form barycentric formula with weights
%   given by the column vector VK to evaluate a polynomial at the points in the
%   column vector X, where the polynomial interpolates the values in the columns
%   of FVALS on the grid of points in [-1, 1] in the column vector XK.
%
% See also CHEBTECH.CHEBPTS, CHEBTECH.BARYWTS, CHEBTECH.FEVAL, CHEBTECH.CLENSHAW.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

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
    warning('CHEBFUN:chebtech:bary:colvec', 'Input should be a column vector.');
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
if ( numel(x) < length(xk) )  % Loop over evaluation points
    % Initialise return value:
    fx = zeros(size(x, 1), m);

    % Loop:
    for j = 1:numel(x),
        xx = vk ./ (x(j) - xk);
        fx(j,:) = (xx.'*fvals) / sum(xx);
    end
else                         % Loop over barycentric nodes
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
for k = find(isnan(fx(:,1)))'
    index = find(x(k) == xk, 1);
    if ( ~isempty(index) )
        fx(k,:) = fvals(index,:);
    end
end

end