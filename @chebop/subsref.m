function varargout = subsref(N, index)
%SUBSREF   Evaluate a CHEBOP or reference its fields
%
% See also: chebop/feval

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
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
        % Evaluate the operator with the input arguments. This method actually
        % does all the tricky set-up for CHEBOP/FEVAL().
        numberOfInputs = nargin(N.op);
        
        if ( numberOfInputs == 1 )
            % Easy case, no X variable involved.
            varargout{1} = feval(N, idx{:});
        else
            if ( numel(idx) > 1 )
                % More than one inputs -- first one must be the X variable
                x = idx{1};
                idx(1) = [];    % Throw away X
            else
                % Create an X variable if it was not passed by the user.
                x = chebfun(@(x) x, N.domain);
            end
            
            % If the second argument is a CHEBMATRIX, but numberOfInputs >
            % 2, we must expand chebmatrix entries out to a cell for {:} to
            % work below.
            
            if ( numberOfInputs > 2 ) && ( isa(idx{1}, 'chebmatrix') )
                args = {};
                for k = 1:numel(idx)
                    if ( isa(idx{k}, 'chebmatrix') )
                        args = [args ; idx{k}.blocks];
                    else
                        args = [args ; idx(k)];
                    end
                end
            else
                args = idx;
            end
            %                 if ( numel(idx) == (numberOfInputs - 1) )
            %                     % Create an X variable if it was not passed by the user.
            %                     x = chebfun(@(x) x, N.domain);
            %                 elseif (  numel(idx) == numberOfInputs )
            %                     % Otherwise, the X variable will be the first element
            %                     x = idx{1};
            %                     idx(1) = []; % Throw away X
            %                 else
            %                     error('CHEBFUN:CHEBOP:subsref', ...
            %                         'Unexpected number of inputs');
            %                 end
            varargout{1} = feval(N, x, args{:});
        end
    otherwise
        
        error('CHEBOP:subsref:indexType',...
            ['Unexpected index.type of ' index(1).type]);
        
end