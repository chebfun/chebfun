function [edge, vscale] = detectedge(op, domain, hscale, vscale) %, der)
%EDGEDETECT   Edge detection.
%
% EDGE = DETECTEDGE(F, A, B, HSCALE, VSCALE, DER) detects a blowup in first,
% second, third, or fourth derivatives of F in [A,B]. HSCALE is the horizontal
% scale and VSCALE is the vertical scale. If no edge is detected, EDGE = 0 is
% returned.
%
% DER is optional and is the derivative of a map (a function handle). It is used
% in the unbounded domain case. If it is not provided, the identity map is
% assumed.

% Call the main routine:
[edge, vscale] = detectedgeMain(op, domain, hscale, vscale, @(x)0*x+1);

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

function [edge, vscale] = detectedgeMain(op, domain, hscale, vscale, der)

% checkblowup = false;

a = domain(1);
b = domain(2);

% Assume no edge is found
edge = [];

nder = 4; % Number of derivatives to be tested
N = 15;   % grid size for finite difference computations in loop.

% Compute norm_inf of first nder derivatives. FD grid size is 50.
[na, nb, maxd] = findMaxDer(op, a, b, nder, 50, der);

% Keep track of endpoints:
ends = [na(nder), nb(nder)];

% Main loop:
while ( maxd(nder) ~= inf && ~isnan(maxd(nder)) &&  diff(ends) > eps*hscale )
    
    % Keep track of previous max derivatives:
    maxd1 = maxd(1:nder);           
    
    % Compute maximum derivatives on interval:
    [na, nb, maxd] = findMaxDer(op, ends(1), ends(2), nder, N, der);  
    
    % Find proper nder: (i.e., how many derivatives to test)
    nder = find( (maxd > (5.5-(1:nder)').*maxd1 ) & ...
        (maxd > 10*vscale./hscale.^(1:nder)') , 1, 'first' );      
    
    if ( isempty(nder) )
        % Derivatives are not growing; return edge = [].
        return                                          
    elseif ( nder == 1 && diff(ends) < 1e-3*hscale )
        % Blow up in first derivative; use findjump().
        edge = findJump(op, ends(1) ,ends(2), hscale, vscale, der);  
        return
    end
    
    % Edge is within this interval. Now repeat process.
    ends = [na(nder),  nb(nder)];
    
%     % Blowup mode? [TODO]: Re-enable blowup mode.
%     if ( checkblowup && abs(op(mean(ends))) > 100*vscale )
%         nedge = findBlowup(f, ends(1), ends(2) , vscale, hscale);
%         if ( isempty(nedge) )
%             checkblowup = false;
%         else
%             edge = nedge;
%             return
%         end
%     end  
    
end

edge = mean(ends);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edge = findJump(op, a, b, hscale, vscale, der)
% Detect blowup in first the derivative and use bisection to locate the edge.

% Assume no edge has been found:
edge = [];                                  

% Compute values at ends:
y = op([a ; b]);
ya = y(1,:); 
yb = y(2,:); 
% yb = op(b);       

% Estimate max abs of the derivative:
maxd = abs( ya - yb ) / ( (b - a).*der(b + a)/2);                   

% If derivative is very small, this is probably a false edge.
if ( maxd < 1e-5 * vscale/hscale )
    return
end

% Keep track how many times derivative stopped growing:
cont = 0;                                   

% Estimate edge location:
e1 = (b + a)/2;                               

% Force loop
e0 = e1 + 1;                                  

% Main loop: (Note that maxd = inf whenever dx < realmin)
while ( ( cont < 2 || maxd == inf ) && ( e0 ~= e1 ) )
    % Evaluate OP at c, the center of the interval [a,b]:
    c = (a+b)/2; 
    yc = op(c);                 
    
    % Find the undivided difference on each side of interval
    dyl = abs(yc-ya) / der((a+c)/2);
    dyr = abs(yb-yc) / der((b+c)/2);          
    
    % Keep track of maximum value:
    maxd1 = maxd;                           
    
    if ( dyl > dyr )
        % Blow-up seems to be in [a,c]. Bisect:
        b = c; 
        yb = yc;                      
        % Update maxd:
        maxd = dyl/(b-a);                       
    else
        % Blow-up seems to be in [c,b]. Bisect:
        a = c; 
        ya = yc;        
        % Update maxd:
        maxd = dyr/(b-a);
    end
    
    % Update edge location:
    e0 = e1;
    e1 = (a+b)/2;
    
    % Test must fail twice before breaking the loop:
    if ( maxd < maxd1*(1.5) )
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [na, nb, maxd] = findMaxDer(op, a, b, nder, N, der)
% Compute the norm_inf of derivatives 1:nder of f.

% Initial setup:
maxd = zeros(nder, 1);
na = a*ones(nder, 1);
nb = b*ones(nder, 1);

% Generate FD grid points and values
dx = (b-a)/(N-1);
x = [a + (0:N-2)*dx, b].';
dy = op(x);

% Main loop (through derivatives), undivided differences
for j = 1:nder
    dy = diff(dy);
    x = ( x(1:end-1) + x(2:end) )/2;
    [maxd(j), ind] = max(abs(dy./der(x)));
    if ( ind > 1)
        na(j) = x(ind-1);
    end
    if ( ind < length(x) - 1)
        nb(j) = x(ind+1);
    end
end
if ( dx^nder <= eps(0) )
    % Avoid divisions by zero!
    maxd = inf + maxd;
else
    % Get norm_inf of derivatives.
    maxd = maxd./dx.^(1:nder)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function edge = findBlowup(f, a, b , vs, hs)
% % Detects blowup location in function values.
%         
% % Evaluate op at the ends of [a, b]:
% x = [a ; b];
% y = abs(f(x));
% ya = y(1);
% yb = y(2); 
% 
% while ( b-a > 1e7*hs )
%     x = linspace(a, b, 50).';
%     yy = abs(f(x(2:end-1)));
%     y = [ya ; yy(:) ; yb];
%     [maxy, ind] = max(abs(y));
%     if ( ind == 1 )
%         b = x(3);     yb = y(3);
%     elseif ( ind == 50 )
%         a = x(48);    ya = y(48);
%     else
%         a = x(ind-1); ya = y(ind-1); % [TODO]: Check this
%         b = x(ind+1); yb = y(ind+1); % [TODO]: Check this
%     end
% end
% 
% while (b - a > 50*eps(a) )
%     x = linspace(a, b, 10).';
%     yy = abs(f(x(2:end-1)));
%     y = [ya; yy(:); yb];
%     [maxy, ind] = max(abs(y));
%     if ( ind == 1 )
%         b = x(3);     yb = y(3);
%     elseif (ind == 10)
%         a = x(8);     ya = y(8);
%     else
%         a = x(ind-1); ya = y(ind-1); % [TODO]: Check this
%         b = x(ind+1); yb = y(ind+1); % [TODO]: Check this
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
