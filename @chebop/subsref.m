function varargout = subsref(f,index)
%SUBSREF   Evaluate a CHEBOP or reference its fields

% Copyright 2011 by The University of Oxford and The Chebfun Developers. 
% See http://www.maths.ox.ac.uk/chebfun/ for Chebfun information.

idx = index(1).subs;
switch index(1).type
    case '.'
        varargout = { get(f,idx) };
        if length(index)>1
          fun = @(v) subsref(v,index(2:end));
          varargout = cellfun(fun,varargout,'uniform',false);
        end
    case '()'
        varargout{1} = feval(f,idx{:});    
    otherwise
        error('CHEBOP:subsref:indexType',...
            ['Unexpected index.type of ' index(1).type]);
end