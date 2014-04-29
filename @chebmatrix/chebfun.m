function f = chebfun(A)
%CHEBFUN Convert a CHEBMATRIX to an array-valued CHEBFUN, if possible. 
%   F = CHEBFUN(A) converts a chebmatrix A to an array-valued chebfun F, if
%   A is a row or column chebmatrix whose entries are all chebfuns or
%   scalars.
%
%   See also CHEBMATRIX.

%  Copyright 2014 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

if ( ~any( size(A)==1 ) )
    error('CHEBFUN:CHEBMATRIX:chebfun', ...
        'Only a row or column chebmatrix can be converted to chebfun.')
end

cls = blockClasses(A);
isFun = strcmp('chebfun',cls);
isVal = strcmp('double',cls);

% Step through the blocks and concatenate horizontally. 
f = [];
for n = 1:length(A.blocks)
    if ( isFun(n) )
        f = [f, A.blocks{n}];
    elseif ( isVal(n) )
        f = [f, chebfun(A.blocks{n}, A.domain)];
    else
        error('Cannot convert element %i to a chebfun.',n)
    end
end

end
