function measure = measure(f, a, b)
%MEASURE    Measure of a CHEBFUN F on an interval.
%   MEASURE(F, A, B) computes the number F^-1([a,b]) i.e., the measure of
%   the set which is mapped to values between A and B under the mapping F.
%   MEASURE(F, [A, B]) is an equivalent syntax.

% Copyright 2016 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

% TODO: This method needs a test.

if ( nargin == 2 )
    d = a;
    a = d(1);
    b = d(2);
else
    d = [a, b];
end

% d must be a 2X1 or 1X2 vector:
if ( ~isvector(d) || length(d) > 2)
    error('CHEBFUN:CHEBFUN:measure:badInput', ...
        'Invalid second input argument.');
end

% Function f must be real valued:
if ( ~isreal(f) )
    error('CHEBFUN:CHEBFUN:measure:complex', ...
        'Complex valued functions can not be measured in Chebfun.');
end

% If f is not continuous, results may be wrong.
if ( ~iscont(f) )
    warning('CHEBFUN:CHEBFUN:measure:discontinuous', ...
        'The function is not continuous, results may be inaccurate.');
end

% a should not be greater than b:
if ( (a > b) || (a == inf) || (b == -inf) )
    error('CHEBFUN:CHEBFUN:measure:badInterval', ...
        'a should not be greater than b.' );
end

% What to do if a == b
% if ( a == b ) end

if ( (a == -inf) && (b == inf) )
    measure = f.ends(end) - f.ends(1);
    return
end

% If a = -inf and b is bounded:
if ( (a == -inf) || (a < min(f)) )
    % A continuous function can not hit a:
    ra = [];
else
    % Compute the points where the function hits a:
    ra = roots(f-a);
end

% If b == inf and a is bounded:
if ( (b == inf) || (b > max(f)) )
    % A continuous function can not hit b:
    rb = [];
else
    % Compute the points where the function hits b:
    rb = roots(f-b);
end

% Sort and 'typify' these points into a single array.
[r, t] = sortAndMerge(ra, rb);

% Compute the derivative of f:
fp = diff(f);

% Compute the measure:
measure = 0;
for i = 1:length(r)-1
    if ( t(i) == 'a' )
        if ( t(i+1) == 'b' )
            % Function hits a then b
            measure = measure + r(i+1) - r(i);
        else
            % Two consecutive roots hitting a, make sure the function is
            % inside a and b:
            if ( (feval(fp, r(i)) >= 0) && (feval(fp, r(i+1)) <= 0) )
                measure = measure + r(i+1) - r(i);
            end
        end
    else
        if ( t(i+1) == 'a' )
            % Function hits b then a
            measure = measure + r(i+1) - r(i);
        else
            % Two consecutive roots hitting b, make sure the function is
            % inside a and b:
            if ( (feval(fp, r(i)) <= 0) && (feval(fp, r(i+1)) >= 0) )
                measure = measure + r(i+1) - r(i);
            end
        end
    end
end

% Search for measure at the boundary points.
ends = [f.domain(1), f.domain(end)];
fLeft = feval(f, ends(1));
fRight = feval(f, ends(end));

if ( ~isempty(r) )
    % Check if there is any measure from the left boundary to the first root.
    if ( (fLeft >= a) && (fLeft <= b) )
        if ( t(1) == 'a' )
            % Function hits a, check derivative:
            if ( feval(fp, r(1)) <= 0 )
                measure = measure + r(1) - ends(1);
            end
        else
            % Function hits b, check derivative:
            if ( feval(fp, r(1)) >= 0 )
                measure = measure + r(1) - ends(1);
            end
        end
    end
    
    % Check if there is any measure from the last root to the right boundary.
    if ( (fRight >= a) && (fRight <= b) )
        if ( t(end) == 'a' )
            % Function hits a, check derivative:
            if ( feval(fp, r(end)) >= 0 )
                measure = measure + ends(end) - r(end);
            end
        else
            % Function hits b, check derivative:
            if ( feval(fp, r(end)) <= 0 )
                measure = measure + ends(end) - r(end);
            end
        end
    end
else
    % Function does not hit a or b.
    if ( (fLeft >= a) && (fLeft <= b) )
        if ( (fRight >= a) && (fRight <= b) )
            % In this case measure is the length of the whole interval.
            measure = ends(end) - ends(1);
        else
            % This means the function did hit a or b and but not detected.
            error('CHEBFUN:CHEBFUN:measure:noRoot', 'Root not detected');
        end
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r, t] = sortAndMerge(ra, rb)
%SORTANDMERGE constructs a single sorted array obtained from ra and rb.
%   SORTANDMERGE(RA, RB) merges RA and RB and sorts the merged list. Vector
%   T is of the same length as the merged array such that T(i) is either
%   'a' or 'b' indicating the original array from which the ith member of
%   the merged array is coming from. It is assumed that RA and RB have no
%   elements in common.

m = 1; n = 1; i = 1;
lenA = length(ra);
lenB = length(rb);

r = zeros(lenA + lenB, 1);
t = zeros(lenA + lenB, 1);

while ( (m <= lenA) || (n <= lenB) )
    % If one of the arrays is exausted, copy the second one into output:
    if ( m > lenA )
        while ( n <= lenB )
            r(i) = rb(n);
            t(i) = 'b';
            n = n + 1;
            i = i + 1;
        end
        break
    end
    
    if ( n > lenB )
        while ( m <= lenA )
            r(i) = ra(m);
            t(i) = 'a';
            m = m + 1;
            i = i + 1;
        end
        break
    end
    
    % None of the arrays is exausted;
    if ( ra(m) < rb(n) )
        r(i) = ra(m);
        t(i) = 'a';
        i = i + 1;
        m = m + 1;
    else
        r(i) = rb(n);
        t(i) = 'b';
        i = i + 1;
        n = n + 1;
    end
end

% Sanity check:
if ( length(r) ~= lenA + lenB )
    % Something has gone wrong.
    error('CHEBFUN:CHEBFUN:measure:sortAndMerge:duplicateRoots', ...
        'Roots likely to have duplication.');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isc = iscont(F)
%ISCONT   Continuity test for chebfuns.
%   ISCONT(F) Returns logical ture if the CHEBFUN F is continuous on the closed
%   interval defined by the domain of F and logical 0 otherwise.

for k = 1:numel(F)
    isc = contCheck(F(k));
    if ( ~isc )
        return
    end
end

end

function isc = contCheck(f)
% Check continuity of a CHEBFUN f.
isc = true;

% Delta functions.
pref = chebfunpref();
deltaTol = pref.deltaPrefs.deltaTol;
out = get(f, 'deltas');
if ( ~isempty(out) )
    deltas = out(2:end, :);
    if ( any(any(abs(deltas(:)) > deltaTol)) )
        % Function with non-trivial delta functions is never continuous.
        isc = false;
        return
    end
end

% Exponents.
if ( ~isfinite(f) )
    % Function with blowup is never continuous.
    isc = false;
    return
end

% Smooth CHEBFUN.
if ( numel(f.funs) == 1 )
    % Function with a single fun is always continuous.
    isc = true;    
    return
end

% Function has more than one FUN:
dom = f.domain;
valTol = 100*max(vscale(f)*eps);
for i = 2:length(dom)-1        
    fLeft = feval(f, dom(i), 'left');
    fRight = feval(f, dom(i), 'right');
    fMiddle = f.pointValues(i, 1);
    if ( any(abs([fLeft-fRight, fLeft-fMiddle, fRight-fMiddle]) > valTol) )
        isc = false;
        return
    end        
end    

end
