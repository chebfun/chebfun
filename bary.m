function fx = bary(x, fvals, xk, vk)
%BARY   Barycentric interpolation formula.
%   BARY(X, FVALS, XK, VK) uses the 2nd form barycentric formula with weights VK
%   to evaluate an interpolant of the data {XK, FVALS(:,k)} at the points X.
%   Note that XK and VK should be column vectors, and FVALS, XK, and VK should
%   have the same length.
%
%   BARY(X, FVALS) assumes XK are the 2nd-kind Chebyshev points and VK are the
%   corresponding barycentric weights.
%
%   If size(FVALS, 2) > 1 then BARY(X, FVALS) returns values in the form
%   [F_1(X), F_2(X), ...], where size(F_k(X)) = size(X) for each k.
%
%   Example:
%     x = chebpts(181);
%     f = 1./( 1 + 25*x.^2 );
%     xx = linspace(-1, 1, 1000);
%     [xx, yy] = meshgrid(xx, xx);
%     ff = bary(xx + 1i*yy, f);
%     h = surf(xx, yy, 0*xx, angle(-ff));
%     set(h, 'edgealpha', 0)
%     view(0, 90), shg
%     colormap(hsv)
%
% See also CHEBTECH.CLENSHAW.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
[n, m] = size(fvals);
sizex = size(x);
ndimsx = ndims(x);

if ( (m > 1) && (ndimsx > 2) )
    error('CHEBFUN:bary:evalArrayAtNDArray', ...
        ['BARY does not support evaluation of vectors of polynomials at ' ...
         'inputs with more than two dimensions.']);
end

% Default to Chebyshev nodes and barycentric weights:
if ( nargin < 3 )
    xk = chebtech2.chebpts(n);
end

if ( nargin < 4 )
    vk = chebtech2.barywts(n);
end

if ( ~all(sizex) )
    fx = x;
    return
end

% Check that input is a column vector:
if ( (ndimsx > 2) || (sizex(2) > 1) )
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

% Reshape the output if possible:
if ( (m == 1) && ( (ndimsx > 2) || (sizex(2) > 1) ) )
    fx = reshape(fx, sizex);
elseif ( (m > 1) && ( (ndimsx == 2) || (sizex(2) > 1) ) )
    fx = reshape(fx, sizex(1), m*numel(x)/sizex(1));
end

end
