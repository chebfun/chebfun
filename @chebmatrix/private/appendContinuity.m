function L = appendContinuity(L)
% Append smoothness constraints at domain breakpoints.
d = getVarDiffOrders(L);
dom = domain(L);
if ( max(d) > 0 ) && ( length(dom) > 2 )
    C = domainContinuity(L,max(d)-1);
    Z = linop.evalAt(dom(1),dom)*(linop.zeros(dom));
    zero = chebmatrix({Z});
    for var = 2:length(d)
        zero = [zero,Z];
    end
    for var = 1:length(d)
        B = zero;
        for m = 0:d(var)-1
            for k = 2:length(dom)-1
                B.blocks{var} = C{m+1,k};
                L = L.bc(B,0);
            end
        end
    end
end
end
