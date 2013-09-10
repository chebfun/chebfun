function f = extractColumns(f, colIdx)

% TODO: Document

if ( isempty(f) )
    return
end

if ( ~isnumeric(colIdx) && strcmp(colIdx, ':') )
    return
end

numCols = min(size(f));
if ( max(colIdx) > numCols )
    error('CHEBFUN:subsref:dimensions','indeex exceeds CHEBFUN dimensions.')    
end

if ( numel(colIdx) == numel(numCols) && all(colIdx == 1:numCols) )
    return
end

for k = 1:numel(f.funs)
    f.funs{k} = extractColumns(f.funs{k}, colIdx);
end

f.impulses = f.impulses(:,colIdx,:);

end


