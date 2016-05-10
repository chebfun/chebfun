function L = deriveContinuity(L, dom, makePeriodic)
%DERIVECONTINUITY Continuity conditions in a piecewise domain.
%   L = DERIVECONTINUITY(L) examines the domain of L and the differential
%   orders of the variables in the system, in order to deduce and encode
%   the appropriate continuity conditions for each variable at every
%   breakpoint. The results are stored in the 'continuity' property of L.
%
%   L has a property called hasGivenJumpAt that signals where automatic
%   continuity conditions are NOT to be applied, because they have been
%   given explicitly in the constraints. 
%
%   L = DERIVECONTINUITY(L,DOMAIN) uses the given domain, merged with
%   L.DOMAIN, in order to derive continuity. In this way you can introduce
%   new breakpoints.
%
%   L = DERIVECONTINUITY(L,DOMAIN,TRUE) ignores breakpoints and enforces
%   appropriate continuity at the boundary points, in order to create a
%   constraint of periodicity. This is called by ADDBC, which then moves
%   the conditions to the 'constraint' property.

%  Copyright 2016 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

if ( nargin < 3 )
    makePeriodic = false;
    if ( nargin < 2 )
        dom = [];
    end
end

dom = domain.merge(dom, L.domain);

diffOrd = L.diffOrder;          % order of each block
diffOrd = max(diffOrd, [], 1);  % max order per variable

cont = L.continuity;            % append, don't overwrite

if ( ( nargin < 2 ) || ~makePeriodic )
    % Use the interior breakpoints for continuity.
    left = dom(2:end-1);
    % Remove any where jumps were given explicitly.
    left = setdiff( left, L.hasGivenJumpsAt ); 
    % Points are the same from both directions. 
    right = left;
else
    % Create periodic conditions, using only endpoints.
    left = dom(end);
    right = dom(1);
end    

if ( ( max(diffOrd) > 0 ) && ( ~isempty(left) ) )
    
    % These give all possible continuity statements up to the maximum
    % order.
    C = domainContinuity(dom,max(diffOrd)-1, left, right);
    
    % Each function variable gets a zero functional block; each scalar variable
    % gets a scalar zero.
    z = functionalBlock.zero(dom);
    Z = {}; %#ok<*AGROW>
    for var = 1:length(diffOrd)
        if ( isnan(diffOrd(var)) || diffOrd(var) == 0 )  % scalar variable
            Z = [Z, 0]; 
        else   % function variable
            Z = [Z, z];
        end
    end
    
    for var = 1:length(diffOrd)
        % Skip if this is a scalar variable; it plays no role in continuity.
        if ( isnan(diffOrd(var)) || diffOrd(var) == 0 )
            continue
        end
        
        B = Z;
        for m = 0:diffOrd(var)-1    % up to this variable's diff. order
            for k = 1:length(left)  % for each point
                B(var) = C{m+1, k}; % right/left difference in mth deriv.
                cont = cont.append(B, 0);
            end
        end
    end
end

L.continuity = cont;

end

function C = domainContinuity(dom, maxorder,left, right)
% Returns expressions of continuity at the breakpoints of the domain of L.
%   C{m,k} has the (m-1)th-order derivative at breakpoint k

% Initialize output cell array
numIntervals = length(left);
C = cell(maxorder+1, numIntervals);

% Construct all evaluation difference functionals, and start by storing those in
% the output cell, C.
for k = 1:numIntervals
    C{1, k} = functionalBlock.feval(left(k), dom, -1) - ...
              functionalBlock.feval(right(k), dom, +1);
end

% Now multiply the top row of C by increasingly high differentation operators:
for m = 1:maxorder
    Dm = operatorBlock.diff(dom, m);
    Dm_rep = repmat({Dm}, 1, numIntervals);
    C(m+1, :) = cellfun(@mtimes, C(1,:), Dm_rep, 'uniformOutput', false);
end 

end
