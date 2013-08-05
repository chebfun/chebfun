function [edge, vscale] = detectEdge(op, domain, vscale, hscale, derHandle)
%DETECTEDGE   Edge detection.
%   EDGE = DETECTEDGE(F, DOMAIN, HSCALE, VSCALE) detects a blowup in first,
%   second, third, or fourth derivatives of F in [A,B]. HSCALE is the horizontal
%   scale and VSCALE is the vertical scale. If no edge is detected, EDGE = 0 is
%   returned.

% Copyright 2013 by The University of Oxford and The Chebfun Developers.
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

%[TODO]: This code will need to be revisited once unbounded domains and blowup
%        are supported again.

%   DERHANDLE is optional and is the derivative of a map (a function handle). It
%   is used in the unbounded domain case. If it is not provided, the identity
%   map is assumed.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PSEUDOCODE. Take from Pachón, Platte and Trefethen, "Piecewise smooth
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
% Parse the inputs:
if ( nargin < 3 )
    hscale = norm(domain, inf);
    if ( isinf(hscale) )
        hscale = 1;
    end
end
if ( nargin < 4 )
    vscale = 0;
else
    vscale = max(vscale);
end
% [TODO]: This may be required when we have unbounded maps again.
if ( nargin < 5 )
    derHandle = @(x) 0*x + 1;
end

% Call the main routine:
[edge, vscale] = detectedgeMain(op, domain, hscale, vscale, derHandle);

% Tidy the results:
% If we didn't detect an edge, then bisect:
if ( isempty(edge) )
    edge = mean(domain);
end
htol = 1e-14*hscale;
% If the edge as at the end of the domain, move it in by 1%:
if ( abs(domain(1) - edge) <= htol )
    edge = domain(1) + diff(domain)/100;
elseif ( abs(domain(2) - edge) <= htol )
    edge = domain(2) - diff(domain)/100;
end

end

function [edge, vscale] = detectedgeMain(op, domain, hscale, vscale, derHandle)

% checkblowup = false;

a = domain(1);
b = domain(2);

% Assume no edge is found
edge = [];

numTestDers = 4;  % Maximum number of derivatives to be tested.
gridSize1 = 50;   % Grid size for 1st finite difference computations.
gridSize234 = 15; % Grid size for higher derivative computations in loop.

% Compute norm_inf of first numTestDers derivatives.
[new_a, new_b, maxDer] = findMaxDer(op, a, b, numTestDers, gridSize1, ...
    derHandle);

% Keep track of endpoints:
ends = [new_a(numTestDers), new_b(numTestDers)];

% Main loop:
while ( (maxDer(numTestDers) ~= inf) && ~isnan(maxDer(numTestDers)) ...
    &&  (diff(ends) > eps*hscale) )

    % Keep track of previous max derivatives:
    maxDerPrev = maxDer(1:numTestDers);

    % Compute maximum derivatives on interval:
    [new_a, new_b, maxDer] = ...
        findMaxDer(op, ends(1), ends(2), numTestDers, gridSize234, derHandle);

    % Choose how many derivatives to test in this iteration:
    numTestDers = find((maxDer > (5.5 - (1:numTestDers)').*maxDerPrev ) & ...
        (maxDer > 10*vscale./hscale.^(1:numTestDers)'), 1, 'first');

    if ( isempty(numTestDers) )
        % Derivatives are not growing; return edge = [].
        return
    elseif ( (numTestDers == 1) && (diff(ends) < 1e-3*hscale) )
        % Blow up in first derivative; use findjump().
        edge = findJump(op, ends(1) ,ends(2), hscale, vscale, derHandle);
        return
    end

    % Edge is within this interval. Now repeat process.
    ends = [new_a(numTestDers),  new_b(numTestDers)];

%     % Blowup mode? [TODO]: Re-enable blowup mode.
%     if ( checkblowup && abs(op(mean(ends))) > 100*vscale )
%         nedge = findBlowup(f, ends(1), ends(2), gridSize1, gridSize234, ...
%             vscale, hscale);
%         if ( isempty(nedge) )
%             checkblowup = false;
%         else
%             edge = nedge;
%             return
%         end
%     end

end

edge = mean(ends);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edge = findJump(op, a, b, hscale, vscale, derHandle)
% Detect blowup in first the derivative and use bisection to locate the edge.

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
    dyl = abs(yc - ya) / derHandle((a + c)/2);
    dyr = abs(yb - yc) / derHandle((b + c)/2);

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

if ( (e0 - e1) <= eps(e0) )
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
function [na, nb, maxDer] = findMaxDer(op, a, b, numTestDers, gridSize, ...
    derHandle)
% Compute the norm_inf of derivatives 1:nder of f.

% Initial setup:
maxDer = zeros(numTestDers, 1);
na = a*ones(numTestDers, 1);
nb = b*ones(numTestDers, 1);

% Generate FD grid points and values
dx = (b - a)/(gridSize - 1);
x = [a + (0:gridSize-2)*dx, b].';
dy = op(x);

% Main loop (through derivatives), undivided differences:
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

% function edge = findBlowup(f, a, b , gridSize1, gridSize234, vs, hs)
% % Detects blowup location in function values.
%
% % Evaluate op at the ends of [a, b]:
% x = [a ; b];
% y = abs(f(x));
% ya = y(1);
% yb = y(2);
%
% while ( b-a > 1e7*hs )
%     x = linspace(a, b, gridSize1).';
%     yy = abs(f(x(2:end-1)));
%     y = [ya ; yy(:) ; yb];
%     [maxy, ind] = max(abs(y));
%     if ( ind == 1 )
%         b = x(3);     yb = y(3);
%     elseif ( ind == gridSize1 )
%         a = x(gridSize1-2);    ya = y(gridSize1-2);
%     else
%         a = x(ind-1); ya = y(ind-1);
%         b = x(ind+1); yb = y(ind+1);
%     end
% end
%
% while (b - a > 50*eps(a) )
%     x = linspace(a, b, gridSize234).';
%     yy = abs(f(x(2:end-1)));
%     y = [ya; yy(:); yb];
%     [maxy, ind] = max(abs(y));
%     if ( ind == 1 )
%         b = x(3);     yb = y(3);
%     elseif (ind == gridSize234)
%         a = x(gridSize234-1);     ya = y(gridSize234-2);
%     else
%         a = x(ind-1); ya = y(ind-1);
%         b = x(ind+1); yb = y(ind+1);
%     end
% end
%
% while ( b - a >= 4*eps(a) )
%     x = linspace(a, b, 4).';
%     yy = abs(f(x(2:end-1)));
%     y = [ya; yy(:); yb];
%     if ( y(2) > y(3) )
%         b = x(3); yb = y(3);
%     else
%         a = x(2); ya = y(2);
%     end
% end
%
% [maxy, ind] = max(y);
% edge = x(ind);
%
% if ( maxy < 1e5*vs )
%     edge = [];
% end
%
% end
