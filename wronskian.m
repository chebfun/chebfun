function w = wronskian(L, varargin)
%WRONSKIAN   Wronskian of chebfuns.
%   WRONSKIAN(f, g, ..., h) computes the wronskian of the chebfuns.
%
% See also

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

n = length(varargin);
A = zeros(n);
for i = 1:n
    f = varargin{i};
    f.domain = 

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
