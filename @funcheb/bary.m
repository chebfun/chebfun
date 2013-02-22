function fx = bary(x, gvals, xk, vk, kind)
%BARY  Barycentric interpolation on a 2nd-kind Chebyshev grid.
%   BARY(X, GVALS, XK, VK) evaluates the barycentric interpolant (with weights
%   VK) of the interpolant (at the points XK) to the values in the columns of
%   GVALS at the points X.
%    
%   If size(GVALS, 2) > 1 then X should be a column vector. If it is not, a
%   warning is displayed and BARY attempts to return values in the form [G_1(X),
%   G_2(X), ...], where size(G_k(X)) = size(X).
%
%   BARY(X, GVALS, XK, VK, KIND) determines which 'kind' of barycentric
%   interpolant to use, where KIND may be on of 1 or 2. By default the 2nd-kind
%   barycentric formula when evaluating within [-1, 1], and the 1st-kind formula
%   for outside the interval or in the complex plane. (See [1] for details).
%
%   BARY(X, GVALS, KIND) overrides the default behaviour and uses the KIND
%   barycentric formula, where KIND may be either 1 or 2. 
%
%   Example:
%     xcheb = funcheb2.chebpts(14);
%     vcheb = funcheb2.barywts(14);
%     fx = 1./( 1 + 25*xcheb.^2 );
%     xx = linspace(-1, 1, 1000);
%     [xx, yy] = meshgrid(xx, xx);
%     ff = bary(xx + 1i*yy, fx, xcheb, vcheb);
%     h = surf(xx, yy, 0*xx, angle(-ff));
%     set(h, 'edgealpha', 0)
%     view(0,90), shg
%
%   [1] Webb, Trefethen, and Gonnet, "Stability of Barycentric interpolation
%   formulas for extrapolation", SIAM J. Sci. Comput., 2012.
%
% See also FUNCHEB.CHEBPTS, FUNCHEB.BARYWTS, FUNCHEB.FEVAL.

%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

% Parse inputs
[n, m] = size(gvals);
sizex = size(x);
ndimsx = ndims(x);

% If true, possibly use type-1 barycentric formula:
bary1flag = true;  

if ( nargin < 3 || isempty(kind))
    kind = [];
else
    % If we're given a KIND, we don't allow a choice.
    bary1flag = false;
end

% Convert arrays of evaluation points to column vector:
if ( (ndimsx > 2) || (sizex(2) > 1) )
    x = x(:);
end

% The function is a constant!
if ( n == 1 )  
    fx = repmat(gvals, length(x), 1);
    return
end

% The function is NaN!
if ( any(isnan(gvals)) )     
    fx = NaN(sizex);
    return
end

% Call a barycentric formula of type 1 or 2:
if ( kind == 1 )
    fx = bary1(x, gvals, xk, vk);
    
elseif ( kind == 2 )
    fx = bary2(x, gvals, xk, vk);
    
elseif ( bary1flag )
    % Call both! (1 for points outside [-1,1], 2 otherwise).
    mask = imag(x) | x < -1 | x > 1;
    ind1 = find(any(mask,2));
    if ( ~isempty(ind1) )
        fx = NaN(size(x, 1), size(gvals, 2));
        fx(ind1,:) = bary1(x(ind1), gvals, xk, vk);
        ind2 = find(any(isnan(fx), 2));
        fx(ind2,:) = bary2(x(ind2), gvals, xk, vk);
    else
        fx = bary2(x, gvals, xk, vk);
    end
    
end

% Try to clean up NaNs:
for k = find(isnan(fx(:,1)))'
    indx = find(x(k) == xk,1);
    if ( ~isempty(indx) )
        fx(k,:) = gvals(indx,:);
    end
end

% Reshape if possible:
if ( ((ndimsx > 2) || (sizex(2) > 1)) && (m == 1) )
    fx = reshape(fx, sizex);
    
elseif ( ((ndimsx == 2) || (sizex(2) > 1)) && (m > 1))
    fx = reshape(fx, sizex(1), size(gvals, 2)*numel(x)/sizex(1));
    
end

end

%% KIND2

function fx = bary2(x, gvals, xk, ek)
% Evaluate the second-kind barycentric formula. Typically this is the standard
% for evaluating a barycentric interpolant on the interval.
m = size(gvals, 2);

if ( numel(x) < length(xk) )  % Loop over evaluation points
    % Initialise return value:
    fx = zeros(size(x, 1), m);  
    % Loop:
    for j = 1:numel(x),
        xx = ek ./ (x(j) - xk);
        fx(j,:) = (xx.'*gvals) / sum(xx);
    end
else                         % Loop over barycentric nodes
    % Initialise:
    num = zeros(size(x, 1), m); 
    denom = num; 
    % Loop:
    for j = 1:length(xk),
        tmp = (ek(j) ./ (x - xk(j)));
        num = num + bsxfun(@times, tmp, gvals(j,:));
        denom = bsxfun(@plus, denom, tmp);
    end
    fx = num ./ denom;
end

end

%% KIND1

function fx = bary1(x, gvals, xk, vk)
% Evaluate the first-kind barycentric formula. Typically we use this formula for
% evaluating outside the interval [-1, 1]. If the number of nodes is >= 600, we
% compute the log of the nodal polynomial in order to avoid under-/overflow.

n = length(xk);
scale = 2; 
x = scale*x; 
xk = scale*xk;
fx = zeros(size(x, 1), size(gvals, 2));  % Initialise return value

if ( numel(x) < n )     % Loop over evaluation points
    for j = 1:numel(x),
        fx(j,:) = (vk./(x(j) - xk)).' * gvals;
    end
else                    % Loop over interpolation nodes
%     e = ones(1, size(gvals, 2));
    for j = 1:n,
        y = vk(j) ./ (x - xk(j));
        fx = fx + bsxfun(@mtimes, y, gvals(j,:));
    end
end

% Evaluate nodal polynomial ell
if ( n < 600 )
    ell = ones(size(fx, 1), 1);
    if ( numel(x) < n ) % Loop over evaluation points
        for j = 1:numel(x)
            ell(j,:) = prod(x(j) - xk);
        end
    else                % Loop over interpolation nodes
        for j = 1:n
            ell = ell .* (x - xk(j));
        end
    end
else
    ell = zeros(size(x));
    if ( numel(x) < n)  % Loop over evaluation points
        for j = 1:numel(x)
            ell(j) = sum(log(x(j) - xk));
        end
    else                % Loop over interpolation nodes
        for j = 1:n
            ell = ell + log(x - xk(j));
        end
    end
    ell = exp(ell);
    if ( isreal(x) && isreal(gvals) && isreal(xk) && isreal(vk) )
        ell = real(ell);
    end
end

% Combine for the result:
fx = bsxfun(@times, fx, ell) * ( 1/(scale*(1-n)) * (-2/scale)^(n-2) );

end
