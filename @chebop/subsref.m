function varargout = subsref(N, index)
%SUBSREF    Evaluate a CHEBOP or reference its fields.
%     ( )
%   N(X, U) and N(U) are equaivalent to FEVAL(N, X, U) and FEVAL(N, U),
%   respectively.
%
%     .
%   F.PROP returns the property PROP of F as defined by GET(F, 'PROP').
%
%   N.OP(X, U) or N.OP(U) are equivalent to N(X, U) and N(U), respectively.
%   N.LBC(U), N.RBC(U), N.BC(U), and N.INIT(U) function similarly.
%
%     {}
%   N{ ... } is not supported.
%
% See also CHEBOP/FEVAL

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

idx = index(1).subs;
switch index(1).type
    
    case '.'

        % Allow expansion of CHEBMATRIX inputs via a call to FEVAL:
        if ( (length(index) > 1) && ...
                any(strcmp(idx, {'lbc', 'rbc', 'bc', 'op', 'init'})) )
            N.op = N.(idx);
            varargout = {feval(N, index(2).subs{:})};
            return
        end
        
        % Allow access to any of the properties of F:
        varargout = { N.(idx) };
        
        % Recurse on SUBSREF:
        if ( length(index) > 1 )
            fun = @(v) subsref(v, index(2:end));
            varargout = cellfun(fun, varargout, 'uniform', false);
        end
        
    case '()'
        
        % Evaluate the operator. Basicially a wrapper for CHEBOP/FEVAL().
        varargout{1} = feval(N, idx{:});

    otherwise
        
        error('CHEBFUN:CHEBOP:subsref:indexType',...
            ['Unexpected index type encountered: ' index(1).type]);
        
end
