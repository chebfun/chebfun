function u = linsolve(L, f, type)

if ( nargin < 3 )
    type = L.discretizationType;
end
obj = type([]);  % to get a static method

dimVals = floor(2.^(3:14)); 
[rowSize, colSize] = blockSizes(L.operator);
isFunVariable = isinf(colSize(1,:));

for dim = dimVals
    
    % Set up the linear system:
    [A, b, dom] = linSystem(L, f, dim, type);
    
    % Solve:
    uDiscrete = A\b;
    
    % Break discrete solution into chunks representing functions and scalars:
    m = colSize(1,:);
    numint = length(dom) - 1;
    m(isFunVariable) = dim*numint; % replace Inf with discrete size
    u = mat2cell(uDiscrete, m, 1);
    
    uFun = u(isFunVariable);
    uVals = cell2mat(uFun.');
    
    % Take an arbitrary linear combination:
    s = 1 ./ (3*(1:size(uVals, 2)));
    v = uVals*s(:);
    
    % Test happiness:
    isDone = true;
    epsLevel = 0;
    for i = 1:numint
        vint = v( (i-1)*dim + (1:dim) );
        [t1, t2] = obj.testConvergence(vint);
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
    funvals = reshape(u{k}, dim, numint); % Piecewise defined.
    f = obj.makeChebfun(num2cell(funvals, 1)', dom); 
%     u{k} = simplify(f, epsLevel);
    u{k} = f;
end

u = chebmatrix(u);

end
