function varargout = subsref(A, sr)
%SUBSREF   Extract part or property of a CHEBMATRIX.
%   A(I,J) returns the slice (submatrix) of A as with an ordinary matrix. The
%   result is a CHEBMATRIX.
%
%   A{I,J} returns a single block as its native type (LINBLOCK, CHEBFUN,
%   double).
%
%   A.(property) returns a property of the CHEBMATRIX.
%
% See also CHEBMATRIX.SUBSASGN.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

sr1 = sr(1);
switch ( sr1.type )
    
    case '()'
        varargout{1} = chebmatrix(subsref(A.blocks, sr1));
        
    case '{}'
        [varargout{1:nargout}] = subsref(A.blocks, sr1);
        
    otherwise
        varargout{1} = A.(sr(1).subs);
        
end

if ( length(sr) > 1 )
    % Recurse on SUBSREF:
    varargout = cellfun(@(v) subsref(v, sr(2:end)), varargout, 'uniformOut', false);
end

end