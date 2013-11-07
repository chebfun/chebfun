function [u,disc] = linsolve(L,f,discType)

[rowSize, colSize] = blockSizes(L.operator);
isFunVariable = isinf(colSize(1,:));

if isa(discType,'function_handle')
    disc = discType(L);  % create a discretization object
    disc = mergeDomains(disc,L,f);
    dimVals = floor(2.^(3:14)); 
else
    disc = discType;
    dimVals = disc.dimension(1);   % TODO: highly suspect!
end

L.constraint.operator.domain = disc.domain;
if isempty(L.continuity)
     % Apply continuity conditions:
     disc = deriveContinuity(disc);
end

numint = disc.numIntervals;

for dim = dimVals
    
    % TODO: Allow different numbers of points in different subdomains
    disc.dimension = repmat(dim,[1 numint]);

    b = disc.rhs(f);

    % Factor the matrix
    if ( isempty(disc.LUFactors) ) || ( length(disc.LUFactors{1}) ~= length(b) )
        A = disc.matrix();
        [P,Q] = lu(A);
        disc.LUFactors = {P,Q};
    else
        P = disc.LUFactors{1};
        Q = disc.LUFactors{2};
    end
       
    % Solve:
    DiscreteSol = Q \ (P\b);
    
    % Break discrete solution into chunks representing functions and scalars:
    m = colSize(1,:);
    m(isFunVariable) = sum(disc.dimension); % replace Inf with discrete size
    u = mat2cell(DiscreteSol, m, 1);
    
    uFun = u(isFunVariable);
    uVals = cell2mat(uFun.');
    
    % Take an arbitrary linear combination:
    s = 1 ./ (3*(1:size(uVals, 2)));
    vVals = uVals*s(:);
    v = mat2cell(vVals,disc.dimension,1);
    
    % Test happiness:
    isDone = true;
    epsLevel = 0;
    for i = 1:numint
        [t1, t2] = disc.testConvergence(v{i});
        isDone = isDone && t1;
        epsLevel = max(epsLevel, t2);
    end
    
    if ( isDone )
        break
    end
end

if ( ~isDone )
    warning('Linear system solution may not have converged.')
end

% The variable u is a cell array with the different components of the solution.
% Because each function component may be piecewise defined, we will loop through
% one by one.
for k = find( isFunVariable )
%    funvals = reshape(u{k}, dim, numint); % Piecewise defined.
    u{k} = disc.toFunction(u{k}); 
%     u{k} = simplify(f, epsLevel);
%    u{k} = f;
end

u = chebmatrix(u);

end
