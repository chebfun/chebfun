function C = domainContinuity(disc,maxorder)
% Returns expressions of continuity conditions at
% the breakpoints of the domain of L.
%   C{m,k} has the (m-1)th-order derivative at breakpoint k

d = disc.domain;
A = linBlock.eye(d);
D = linBlock.diff(d,1);
for m = 0:maxorder
    for k = 2:length(d)-1
        El = linBlock.feval(d(k),d,'-');
        Er = linBlock.feval(d(k),d,'+');
        C{m+1,k} = (El-Er)*A;
    end
    A = D*A;
end
end
