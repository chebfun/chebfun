function fx = trigBary(x, fvals, xk, dom)
%TRIGBARY   Trigonometric barycentric interpolation formula.
%   TRIGBARY(X, FVALS, XK, dom) uses the 2nd form barycentric formula 
%   to evaluate an interpolant of the data {XK, FVALS(:,k)}
%   at the points X. The interpolant is supposed to live on the domain 
%   specified in dom. Note that XK should be column vector, and 
%   FVALS should have the same length.
%
%   TRIGBARY(X, FVALS) assumes XK are equally spaced points in [-pi, pi).
%
%   If size(FVALS, 2) > 1 then TRIGBARY(X, FVALS) returns values in the form
%   [F_1(X), F_2(X), ...], where size(F_k(X)) = size(X) for each k.
%
% REFERENCES:
%   [1] Berrut, Jean-Paul. "Baryzentrische Formeln zur trigonometrischen 
%   Interpolation (I)." Zeitschrift für angewandte Mathematik und Physik 
%   ZAMP 35.1 (1984): 91-105.
% 
%   [2] Henrici, Peter. "Barycentric formulas for interpolating trigonometric 
%   polynomials and their conjugates." Numerische Mathematik 33.2 
%   (1979): 225-234.
%
% See also BARY(), TRIGBARYWEIGHTS()

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse inputs:
[n, m] = size(fvals);
sizex = size(x);
ndimsx = ndims(x);

if ( (m > 1) && (ndimsx > 2) )
    error('CHEBFUN:trigBary:evalArrayAtNDArray', ...
        ['TRIGBARY does not support evaluation of vectors of polynomials at ' ...
         'inputs with more than two dimensions.']);
end

if ( min(size(xk)) > 1 )
    error('CHEBFUN:trigBary:xVector', 'xk must be a vector.');
end

if ( nargin < 4 )
    dom = [-pi, pi];
end

% Default to equispaced nodes in [-pi, pi) and barycentric weights:
if ( nargin < 3 )
    % Default points:
    xk = trigtech.trigpts(n, [-pi, pi]);
end

% Map the points to [-pi, pi]
a = dom(1);
b = dom(2);
xk = pi/(b-a)*(2*xk-a-b);
x  = pi/(b-a)*(2*x-a-b);

xkMax = max(xk);
xkMin = min(xk);
if ( (xkMax > pi + 10*eps) || (xkMin < -pi - 10*eps) )
     error('CHEBFUN:trigbary:invalidNodes', 'nodes XK must lie within the domain');
end

% Remove a periodic end-point and interpolate the average:
if ( norm([xk(1), xk(end)] - [-pi, pi]) < 2*pi*eps )
    xk(end) = [];
    fvals(1, :) = (fvals(1, :) + fvals(end, :))/2;
    fvals(end, :) = [];
end

% Evaluate the weights:
vk = trigBaryWeights(xk);

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

% Choose the appropriate function based on the length of the values to be
% interpolated:
if ( rem(n, 2) == 0 )
    s = sum(xk);
    if ( abs(rem(s, pi)) < 4*pi*eps )
        c = 0;
    else        
        c = cot(s/2);
    end    
    ctsc = @(x) cot(x) + c;
else
    ctsc = @(x) csc(x);
end

% The main loop:
% Initialise:
num = zeros(size(x, 1), m);
denom = num;

% Loop:
for j = 1:length(xk),
    tmp = vk(j) * ctsc((x - xk(j))/2);
    num = num + bsxfun(@times, tmp, fvals(j,:));
    denom = bsxfun(@plus, denom, tmp);
end
fx = num ./ denom;

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
