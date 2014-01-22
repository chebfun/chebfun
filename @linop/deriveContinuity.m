function L = deriveContinuity(L,makePeriodic)

% Find automatic smoothness constraints at domain breakpoints.
d = L.diffOrder;
d = max(d,[],1);
dom = L.domain;

cont = L.continuity;

if ( ( nargin < 2 ) || ~makePeriodic )
    % Use the interior breakpoints.
    left = dom(2:end-1);
    right = left;
else
    % Create periodic conditions.
    left = dom(end);
    right = dom(1);
end
    

if ( max(d) > 0 ) && ( ~isempty(left) )
    
    C = domainContinuity(dom,max(d)-1,left,right);
    
    % Each function variable gets a zero functional block; each scalar variable
    % gets a scalar zero.
    z = functionalBlock.zeros(dom);
    Z = {};
    for var = 1:length(d)
        if ( isnan(d(var)) || d(var) == 0 ) % scalar
            Z = [Z,0];
        else
            Z = [Z,z];
        end
    end
    %    Z = chebmatrix(Z);
    
    for var = 1:length(d)
        % Skip if this is a scalar variable; it plays no role in continuity.
        if ( isnan(d(var)) || d(var) == 0 )
            continue
        end
        B = Z;
        for m = 0:d(var)-1
            for k = 1:length(left)
                B(var) = C{m+1,k};
                cont = cont.append(B,0);
            end
        end
    end
end

L.continuity = cont;

end


function C = domainContinuity(dom,maxorder,left,right)
% Returns expressions of continuity conditions at
% the breakpoints of the domain of L.
%   C{m,k} has the (m-1)th-order derivative at breakpoint k

A = operatorBlock.eye(dom);
D = operatorBlock.diff(dom,1);
for m = 0:maxorder
    for k = 1:length(left)
        El = functionalBlock.feval(left(k),dom,-1);
        Er = functionalBlock.feval(right(k),dom,+1);
        C{m+1,k} = (El-Er)*A;
    end
    A = D*A;
end

end
