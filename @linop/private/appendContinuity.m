function L = appendContinuity(L)

% Append smoothness constraints at domain breakpoints.
d = L.blockDiffOrders;
d = max(d,[],1);
dom = L.domain;

if ( max(d) > 0 ) && ( length(dom) > 2 )

    C = domainContinuity(L,max(d)-1);
    
    % Each function variable gets a zero functional block; each scalar variable
    % gets a scalar zero. 
    z = linBlock.zero(dom);
    Z = {};
    for var = 1:length(d)
        if isnan(d(var))  % scalar
            Z = [Z,0];
        else
            Z = [Z,z];
        end
    end
%    Z = chebmatrix(Z);
    
    for var = 1:length(d)
        % Skip if this is a scalar variable; it plays no role in continuity.
        if isnan(d(var))
            continue
        end
        B = Z;
        for m = 0:d(var)-1
            for k = 2:length(dom)-1
                B.blocks{var} = C{m+1,k};
                L = L.addbc(B,0);
            end
        end
    end
end

end
