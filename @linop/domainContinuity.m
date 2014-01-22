function C = domainContinuity(L,maxorder,doPeriodic)
% Returns expressions of continuity conditions at
% the breakpoints of the domain of L.
%   C{m,k} has the (m-1)th-order derivative at breakpoint k

d = L.domain;

if ( nargin < 3 )
    left = d(2:end);
    right = left;
else
    left = d(end);
    right = d(1);
end

A = operatorBlock.eye(d);
D = operatorBlock.diff(d,1);
for m = 0:maxorder
    for k = 1:length(left)
        El = functionalBlock.feval(left(k),d,'-');
        Er = functionalBlock.feval(right(k),d,'+');
        C{m+1,k} = (El-Er)*A;
    end
    A = D*A;
end

end
