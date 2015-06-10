function edge = detectEdge(f, op, hscale, vscale, pref)
%DETECTEDGE   Edge detection.
%   EDGE = DETECTEDGE(F, OP, HSCALE, VSCALE, BLOWUP) detects a blowup in the 
%   first, second, third, or fourth derivatives of OP in the domain of the FUN 
%   F. HSCALE is the horizontal scale and VSCALE is the vertical scale (note 
%   that both are required). If no edge is detected, EDGE is set to the midpoint
%   of DOMAIN.

% Copyright 2015 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSEUDOCODE. Taken from Pachon, Platte and Trefethen, "Piecewise smooth
% chebfuns" (IMA J. Numer. Anal., 2010)

%  function edge = detectedge(f,a,c) % Find singularity of f in [a; c]
%  edge = NaN
%  On a 50-point equispaced grid in [a;c] compute estimates of |f'|,..,|f''''|.
%  Set b = the gridpoint associated with the maximum estimate of |f''''|.
%  Set dmax = 4.
%  while the current interval [a;c] is larger than machine precision
%    Set a and c to the gridpoints left and right of b, respectively.
%    Refine 7-fold to a 15-point grid in [a;c], and find the gridpoints 
%     associated there with the maximum estimates of |f'|,..,|f^(dmax)|
%    if refinement has not increased the estimates by a factor of 1:2 or more
%      return (i.e., no edge has been detected)
%    end if
%    Set dmax = order of lowest derivative that has increased by such a factor.
%    Set b = gridpoint associated with maximum value of |f^(dmax)|
%    if dmax = 1
%      b = findjump(f,a,c)
%      return
%    end if
%  end while
%  edge = b

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the domain:
dom = f.domain;

% We do not allow an edge at the edge of unbounded domains:
allowEdgeAtBdy = isfinite(dom);

% TODO: Adjusting the operator to include exponents as is done below is a clear
% violation of encapsulation. We could get around this by having a
% 'modifyOperator()' method (better ideas for a name are welcome) which does
% this adjustment (including the inclusion of the map in the unbounded case) at
% the lower levels. Ideally this code should live at the SMOOTHFUN level?
exps = get(f, 'exponents');

% Locate the edges/splitting locations:
if ( all( isfinite( dom ) ) )  
    % bounded domain
    forHandle = @(x) x;
    derHandle = @(x) 0*x + 1;
    if ( any(exps) )
        % Compensating for exponents:
        op = @(x) op(x) ./ ((x - dom(1)).^exps(1) .* (dom(2) - x).^exps(2));
    end  
    
else
    % unbounded domain
    forHandle = f.mapping.For;
    derHandle = f.mapping.Der;
    dom = [-1+eps, 1-eps];
    op = @(x) op(forHandle(x));
    if ( any(exps) )
        % nonzero exponents
        op = @(x) op(x) .* ((x + 1).^exps(1) .* (1 - x).^exps(2));
    end
    
end

% Take the maximum of the vscales if a vector is given:
vscale = max(vscale);

% Call the main routine:
edge = detectedgeMain(op, dom, vscale, hscale, derHandle, pref);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Tidy the results  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

htol = 1e-14*hscale;
% If the edge is at the end of the domain, move it in by 1% (unless at +/- inf):
if ( abs(dom(1) - edge) <= htol )
    if ( allowEdgeAtBdy(1) )
        edge = dom(1) + diff(dom)/100;
    else
        edge = []; % Bisect if at infinity.
    end
elseif ( abs(dom(2) - edge) <= htol )
    if ( allowEdgeAtBdy(2) )
        edge = dom(2) - diff(dom)/100;
    else
    	edge = []; % Bisect if at infinity.
    end
end

% If we didn't detect an edge, then bisect:
if ( isempty(edge) )
    edge = mean(dom);
end

% Map the edge back to the physical domain:
edge = forHandle(edge);

end

function edge = detectedgeMain(op, dom, vscale, hscale, derHandle, pref)
%   DERHANDLE is optional and is the derivative of a map (a function handle). It
%   is used in the unbounded domain case. If it is not provided, the identity
%   map is assumed.

if ( nargin < 5 )
    derHandle = @(x) 0*x + 1;
end

% Assume no edge is found:
edge = [];

numTestDers = 4;  % Maximum number of derivatives to be tested.
gridSize1 = 50;   % Grid size for 1st finite difference computations.
gridSize234 = 15; % Grid size for higher derivative computations in loop.

% Compute norm_inf of first numTestDers derivatives.
[new_a, new_b, maxDer] = findMaxDer(op, dom, numTestDers, gridSize1, ...
    derHandle);

% Keep track of endpoints:
ends = [new_a(numTestDers), new_b(numTestDers)];

% Main loop:
checkBlowUp = pref.blowup;
while ( (maxDer(numTestDers) ~= inf) && ~isnan(maxDer(numTestDers)) ...
    &&  (diff(ends) > eps*hscale) )

    % Keep track of previous max derivatives:
    maxDerPrev = maxDer(1:numTestDers);

    % Compute maximum derivatives on interval:
    [new_a, new_b, maxDer] = ...
        findMaxDer(op, ends, numTestDers, gridSize234, derHandle);
    % Choose how many derivatives to test in this iteration:
    numTestDers = find((maxDer > (5.5 - (1:numTestDers)').*maxDerPrev ) & ...
        (maxDer > 10*vscale./hscale.^(1:numTestDers)'), 1, 'first');

    if ( isempty(numTestDers) )
        % Derivatives are not growing; return edge = [].
        return
    elseif ( (numTestDers == 1) && (diff(ends) < 1e-3*hscale) )
        % Blow up in first derivative; use findjump().
        edge = findJump(op, ends, vscale, hscale, derHandle);
        return
    end

    % Edge is within this interval. Now repeat process.
    ends = [new_a(numTestDers),  new_b(numTestDers)];

    % If the function value turns out to be very large between ENDS, then there 
    % might be a blowing up point. 
    if ( checkBlowUp && ( abs(op(mean(ends))) > 1e2*vscale ) )
        blowUpPoint = findBlowup(op, ends, gridSize1, gridSize234, vscale);
        if ( isempty(blowUpPoint) )
            checkBlowUp = false;
        else
            edge = blowUpPoint;
            return
        end
    end

end

edge = mean(ends);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edge = findJump(op, dom, vscale, hscale, derHandle)
% Detect blowup in first the derivative and use bisection to locate the edge.

a = dom(1);
b = dom(2);

% Assume no edge has been found:
edge = [];

% Compute values at ends:
y = op([a ; b]);
ya = y(1,:);
yb = y(2,:);

% Estimate max abs of the derivative:
maxDer = abs(ya - yb) / ((b - a).*derHandle(b + a)/2);

% If derivative is very small, this is probably a false edge.
if ( maxDer < 1e-5*vscale/hscale )
    return
end

% Keep track how many times derivative stopped growing:
cont = 0;

% Estimate edge location:
e1 = (b + a)/2;

% Force loop
e0 = e1 + 1;

% Main loop: (Note that maxd = inf whenever dx < realmin)
while ( ((cont < 2) || any(maxDer == inf)) && (e0 ~= e1) )
    % Evaluate OP at c, the center of the interval [a,b]:
    c = (a + b)/2;
    yc = op(c);

    % Find the undivided difference on each side of interval
    dyl = max(abs(yc - ya) / derHandle((a + c)/2));
    dyr = max(abs(yb - yc) / derHandle((b + c)/2));

    % Keep track of maximum value:
    maxd1 = maxDer;

    if ( dyl > dyr )
        % Blow-up seems to be in [a,c]. Bisect:
        b = c;
        yb = yc;
        % Update maxd:
        maxDer = dyl/(b - a);
    else
        % Blow-up seems to be in [c,b]. Bisect:
        a = c;
        ya = yc;
        % Update maxd:
        maxDer = dyr/(b - a);
    end

    % Update edge location:
    e0 = e1;
    e1 = (a + b)/2;

    % Test must fail twice before breaking the loop:
    if ( maxDer < maxd1*(1.5) )
        cont = cont + 1;
    end
    
end

if ( (e0 - e1) <= 2*eps(e0) )
    % Look at the floating point at the right:
    yright = op(b + eps(b));
    % If there is a small jump, that is it!
    if ( abs(yright - yb) > eps*100*vscale )
        edge = b;
    else
        edge = a;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [na, nb, maxDer] = findMaxDer(op, dom, numTestDers, gridSize, derHandle)
% Compute the norm_inf of derivatives 1:nder of f.

a = dom(1);
b = dom(2);

% Initial setup:
maxDer = zeros(numTestDers, 1);
na = a*ones(numTestDers, 1);
nb = b*ones(numTestDers, 1);

% Generate FD grid points and values
dx = (b - a)/(gridSize - 1);
x = [a + (0:gridSize-2)*dx, b].';
y = op(x);

% Main loop (through derivatives), undivided differences:
dy = y;
for j = 1:numTestDers
    dy = diff(dy);
    x = (x(1:end-1) + x(2:end))/2;
    dydH = max(abs(bsxfun(@rdivide, dy, derHandle(x))), [], 2);
    [maxDer(j), ind] = max(dydH);
    if ( ind > 1)
        na(j) = x(ind-1);
    end
    if ( ind < length(x) - 1)
        nb(j) = x(ind+1);
    end
end

if ( dx^numTestDers <= eps(0) )
    % Avoid divisions by zero!
    maxDer = inf + maxDer;
else
    % Get norm_inf of derivatives.
    maxDer = maxDer./dx.^(1:numTestDers)';
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function blowUpPoint = findBlowup(op, dom, gridSize1, gridSize234, vscale)
% Detects blowup location in function values.
a = dom(1);
b = dom(2);

% Evaluate OP at the ends of [a, b]:
x = [a ; b];
y = abs(op(x));
ya = y(1);
yb = y(2);

% If the current interval is large:
while ( (b - a) > 1e7*eps(a) )
    [a, b, ya, yb] = zoomIn(op, a, b, ya, yb, gridSize1);
end

% If the interval is now of medium size:
while ( (b - a) > 50*eps(a) )
    [a, b, ya, yb] = zoomIn(op, a, b, ya, yb, gridSize234);
end

% If the interval is small now:
while ( (b - a) >= 4*eps(a) )
    
    % Sampling points:
    x = linspace(a, b, 4).';
    
    % Get the function values at the sampling points: 
    y = abs(op(x(2:end-1)));
    y = [ya; y(:); yb];
    
    % Compare the function values at the two interior points:
    if ( y(2) > y(3) )
        b = x(3); 
        yb = y(3);
    else
        a = x(2); 
        ya = y(2);
    end
end

% Locate the maximum function value, which is the blow-up point:
[maxy, ind] = max(y);
blowUpPoint = x(ind);

% However, if the function at the edge is not that big then it might be a fake
% blow-up point:
if ( maxy < 1e5*vscale )
    blowUpPoint = [];
end

end

function [a, b, ya, yb] = zoomIn(op, a, b, ya, yb, gridSize)
% Narrow down the interval in which the blow-up point is located.

% Sampling points:
x = linspace(a, b, gridSize).';

% Get the function values at the sampling points:
y = abs(op(x(2:end-1)));
y = [ya ; y(:) ; yb];

% Locate the maximum point:
[ignored, ind] = max(abs(y));

% Check where the maximum occurs:
if ( ind == 1 )
    % If the maximum occurs at the left endpoint:
    b = x(3);
    yb = y(3);
elseif ( ind == gridSize )
    % If the maximum occurs at the right endpoint:
    a = x(gridSize-2);
    ya = y(gridSize-2);
else
    % If the maximum occurs at any interior points:
    a = x(ind-1);
    ya = y(ind-1);
    b = x(ind+1);
    yb = y(ind+1);
end

end
