function varargout = subsref(N, index)
%SUBSREF   Evaluate a CHEBOP or reference its fields
%
% See also CHEBOP/FEVAL

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% TODO: Improve documentaiton.

idx = index(1).subs;
switch index(1).type
    
    case '.'
        
        if ( any(strcmp(idx, {'lbc', 'rbc', 'bc'})) )
            fun = N.(idx);
            if ( length(index) > 1)
                varargout = {feval(fun, index(2).subs{:})};
            else
                varargout = {fun};
            end
            return
        end
        
        if ( strcmp(idx, 'op') && numel(index) > 1 && strcmp(index(2).type, '()') )
            varargout{1} = feval(N, index(2).subs{:});
            return
        end
        
        varargout = { get(N,idx) };
        if ( length(index) > 1 )
            fun = @(v) subsref(v, index(2:end));
            varargout = cellfun(fun, varargout, 'uniform', false);
        end
        
    case '()'
        
        % Evaluate the operator. Basicially a wrapper for CHEBOP/FEVAL().
        varargout{1} = feval(N, idx{:});

    otherwise
        
        error('CHEBOP:subsref:indexType',...
            ['Unexpected index.type of ' index(1).type]);
        
end