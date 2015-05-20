function B = subsref(A, sr)
%SUBSREF   Extract part or property of a linop.
%   A(I, J) returns the slice (submatrix) of A as with an ordinary matrix.
%   The result is a CHEBMATRIX.
%
%   A{I, J} returns a single block as its native type (linBlock, chebfun,
%   double).
%
%   A.(property) returns a property of the chebmatrix.
%
% See also CHEBMATRIX.SUBSASGN.

%  Copyright 2015 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org/ for Chebfun information.

% The only reason this method exists is to provide backward compatibility with
% the A(n) realization syntax from version 4, which is referenced in
% Approximation Theory and Approximation Practice. 

if ( length(sr)==1 && strcmp(sr(1).type, '()') && length(sr(1).subs)==1 )
    B = feval(A, sr(1).subs{1});
else
    B = subsref@chebmatrix(A, sr);
end
    
end
