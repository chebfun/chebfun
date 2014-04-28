function varargout = subsref(N, index)
%SUBSREF   Evaluate a CHEBOP or reference its fields

% Copyright 204 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

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
        
        varargout = { get(N,idx) };
        if ( length(index) > 1 )
          fun = @(v) subsref(v, index(2:end));
          varargout = cellfun(fun, varargout, 'uniform', false);
        end
        
    case '()'
        
        % Must expand chebmatrix entries out to a cell for {:} to work below.
        args = {};
        for k = 1:numel(idx)
            if ( isa(idx{k}, 'chebmatrix') )
                idx{k}.blocks
                args = [args ; idx{k}.blocks];
            else
                args = [args ; idx(k)];
            end
        end
        varargout{1} = feval(N, args{:});    
        
    otherwise
        
        error('CHEBOP:subsref:indexType',...
            ['Unexpected index.type of ' index(1).type]);
        
end