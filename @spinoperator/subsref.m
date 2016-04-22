function varargout = subsref(S, index)
%SUBSREF    Evaluate a SPINOPERATOR or reference its fields.
%   S(U) is equaivalent to FEVAL(S, U).
%     
%   S.PROP returns the property PROP of S.
%
% See also SPINOPERATOR/FEVAL.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

idx = index(1).subs;
switch index(1).type
    
    case '.'
        
        % Allow access to any of the properties of S:
        varargout = {S.(idx)};
        
    case '()'
        
        % Evaluate the operator, i.e., wrapper for SPINOPERATOR/FEVAL:
        varargout = {feval(S, idx{:})};
        
    otherwise
        
        error('SPINOPERATOR:subsref:indexType', ...
            ['Unexpected index type encountered: ' index(1).type]);
        
end