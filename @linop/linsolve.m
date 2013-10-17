function u = linsolve(L,f,type)

if ( nargin < 3 )
    type = linBlock.defaultDiscretization;
end

dimVals = floor(2.^(5:.5:14));
[rowSize,colSize] = blockSizes(L.operator);
isFunVariable = isinf(colSize(1,:));

for k = 1:length(dimVals)
    dim = dimVals(k);
    
    [A,b,dom] = linSystem(L,f,dim,type);
    numint = length(dom)-1;
    
    uDiscrete = A\b;
    
    % Break the discrete solution into chunks representing
    % functions and scalars.
    m = colSize(1,:);
    m(isFunVariable) = dim*numint; % replace Inf with discrete size
    u = mat2cell(uDiscrete,m,1);
    
    uFun = u(isFunVariable);
    uVals = cell2mat(uFun');
    
    % Take an arbitrary linear combination.
    s = 1 ./ (3*(1:size(uVals,2)));
    v = uVals*s(:);
    
    obj = type([]);  % to get a static method
    isDone = true;
    epsLevel = 0;
    for i = 1:numint
        vint = v( (i-1)*dim + (1:dim) );
        [t1,t2] = obj.testConvergence(vint);
        isDone = isDone && t1;
        epsLevel = max(epsLevel,t2);
    end
    if isDone
        break
    end
end

if ~isDone
    warning('Linear system solution may not have converged.')
end

if ( isa(obj, 'blockUS') )
    u = full(u{1});
    f = chebtech2({[], flipud(u)});
    f = bndfun(f, dom);
    f = chebfun({f});
    u = chebmatrix({f});
    return
end

% The variable u is a cell array with the different components of
% the solution. Because each function component may be piecewise
% defined, we will loop through one by one.
for k = find( isFunVariable )
    funvals = reshape( u{k}, dim, numint );
    f = chebfun( num2cell(funvals,1)', dom );  % piecewise defined
    u{k} = simplify(f, epsLevel );
end

u = chebmatrix(u);

end
