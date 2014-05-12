function result = parsebc(N, val)

% TODO: Document.

numIn = nargin(N);

if ( isempty(val) )
    result = [];
    
elseif ( isnumeric(val) )
    if ( numIn > 2 )
        % Only allow scalar numerical values to be passed if we are dealing with
        % a scalar problem.
        error('CHEBFUN:CHEBOP:PARSEBC', ...
            'Can only assign scalar BCs to scalar problems');
    else
        result = @(u) u - val;
    end
    
elseif ( isa(val, 'function_handle') )
    % If we are dealing with a scalar problem where the independent variable is
    % not specified in the function handle arguments, allow also passing an
    % input function handle that takes one argument. Otherwise, we request that
    % the number of input to the BC function handle is one less than the number
    % of arguments to the OP part.
    if ( ( (numIn <= 1) && (nargin(val) == 1) ) || ...
            (nargin(val) == (numIn - 1)) )
        result = val;
    else
        error('CHEBFUN:CHEBOP:PARSEBC', ...
            'Number of inputs to BCs do not match operator.');
    end
    
elseif ( iscell(val) ) && ( isnumeric(val{1}) )
    % For backward compatability.
    if ( strcmpi(val{2}, 'neumann') && ( numIn <= 2 ) )
        result = @(u) diff(u) - val{1};
    elseif ( strcmpi(val{2}, 'dirichlet') && ( numIn <= 2 ) )
        result = @(u) u - val{1};
    else
         error('CHEBFUN:CHEBOP:PARSEBC', ...
            'Unable to parse cell input to set BC.');
    end
    
elseif ( strcmpi(val, 'neumann') )
    if ( numIn <= 2 )
        result = @(u) diff(u);
    else
        error('CHEBFUN:CHEBOP:PARSEBC', ...
            'Can only assign scalar BCs to scalar problems.');
    end
    
elseif ( strcmpi(val, 'dirichlet') )
    if ( numIn <= 2 )
        result = @(u) u;
    else
        error('CHEBFUN:CHEBOP:PARSEBC', ...
            'Can only assign scalar BCs to scalar problems.');
    end
    
else
    error('CHEBFUN:CHEBOP:PARSEBC', 'Unsupported format of BCs.')

end

end
