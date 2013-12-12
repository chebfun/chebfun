function disc = deriveContinuity(disc)

% Find automatic smoothness constraints at domain breakpoints.
L = disc.source;
d = L.diffOrder;
d = max(d,[],1);
dom = disc.domain;

cont = linopConstraint();

if ( max(d) > 0 ) && ( length(dom) > 2 )
    
    C = domainContinuity(disc,max(d)-1);
    
    % Each function variable gets a zero functional block; each scalar variable
    % gets a scalar zero.
    z = linBlock.zero(dom);
    Z = {};
    for var = 1:length(d)
        if isnan(d(var)) || d(var) == 0 % scalar
            Z = [Z,0];
        else
            Z = [Z,z];
        end
    end
    %    Z = chebmatrix(Z);
    
    for var = 1:length(d)
        % Skip if this is a scalar variable; it plays no role in continuity.
        if isnan(d(var)) || d(var) == 0
            continue
        end
        B = Z;
        for m = 0:d(var)-1
            for k = 2:length(dom)-1
                B(var) = C{m+1,k};
                cont = cont.append(B,0);
            end
        end
    end
end

disc.source.continuity = cont;

end
