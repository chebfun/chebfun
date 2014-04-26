function varargout = subsref(N, index)
%SUBSREF   Evaluate a CHEBOP or reference its fields

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

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
        if length(index)>1
          fun = @(v) subsref(v,index(2:end));
          varargout = cellfun(fun,varargout,'uniform',false);
        end

        
    case '()'
        varargout{1} = feval(N,idx{:});    
    otherwise
        error('CHEBOP:subsref:indexType',...
            ['Unexpected index.type of ' index(1).type]);
end