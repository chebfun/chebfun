function f = chebfun(A)
%CHEBFUN Convert a row chebmatrix to chebfun, if possible. 

if ( size(A,1)~=1 )
    error('Must use a row chebmatrix for conversion to chebfun.')
end

cls = cellfun(@class,A.blocks,'uniform',false);
isFun = strcmp('chebfun',cls);
isVal = strcmp('double',cls);
f = [];
for n = 1:length(A.blocks)
    if ( isFun(n) )
        f = [f, A.blocks{n}];
    elseif ( isVal(n) )
        f = [f, chebfun(A.blocks{n},A.domain)];
    else
        error('Cannot convert element %i to a chebfun.',n)
    end
end
