function result = parseBC(N, BC)
%PARSEBC  Parse boundary conditions for CHEBOP object.
%   This method is not intended for end users. For information about boundary
%   conditions in chebops, see CHEBOP.
%
%   This method is invoked by the set methods for LBC, RBC, and BC. The types
%   and meanings of the allowed input are described in the documentation for
%   CHEBOP. The result is either empty or a function handle that represents the
%   given condition as needed internally.

% Copyright 2014 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

numIn = nargin(N);

if ( isempty(BC) )
    result = [];
    
elseif ( isnumeric(BC) )
    % This means a Dirichlet condition set to the given value. 
    if ( numIn > 2 )
        % Allow only if we are dealing with a scalar problem.
        error('CHEBFUN:chebop:parsebc:numeric', ...
            'Can only assign scalar BCs to scalar problems');
    else
        result = @(u) u - BC;
    end
    
elseif ( isa(BC, 'function_handle') )
    % If we are dealing with a scalar problem where the independent variable is
    % not specified in the function handle arguments, allow also passing an
    % input function handle that takes one argument. Otherwise, we request that
    % the number of input to the BC function handle is one less than the number
    % of arguments to the OP part.
    if ( ( (numIn <= 1) && (nargin(BC) == 1) ) || (nargin(BC) == (numIn - 1)) )
        result = BC;
    else
        error('CHEBFUN:chebop:parsebc:inputs', ...
            'Number of inputs to BCs do not match operator.');
    end
    
elseif ( strcmpi(BC, 'neumann') )
    % Homogeneous Neumann.
    if ( numIn <= 2 )
        result = @(u) diff(u);
    else
        error('CHEBFUN:chebop:parsebc:neuman', ...
            'Can only assign scalar BCs to scalar problems.');
    end
    
elseif ( strcmpi(BC, 'dirichlet') )
    % Homogeneous Dirichlet. 
    if ( numIn <= 2 )
        result = @(u) u;
    else
        error('CHEBFUN:chebop:parsebc:dirichlet', ...
            'Can only assign scalar BCs to scalar problems.');
    end
    
elseif ( iscell(BC) ) && ( length(BC) == 2 ) && ( ischar(BC{2}) ) && ...
        ( isnumeric(BC{1}) )
    % A cell may have a numerical value followed by a keyword of 'dirichlet' or
    % 'neumann'.
    
    % This behavior is retained only for backward compatability and only for
    % problems with one variable. 
    warning('CHEBFUN:chebop:parsebc:keywordbc',...
        ['Keyword/value specifications of boundary conditions are ', ...
         'deprecated and may be removed in future versions of Chebfun.'])
            
    if ( strcmpi(BC{2}, 'neumann') && ( numIn <= 2 ) )
        result = @(u) diff(u) - BC{1};
    elseif ( strcmpi(BC{2}, 'dirichlet') && ( numIn <= 2 ) )
        result = @(u) u - BC{1};
    else
         error('CHEBFUN:chebop:parsebc:cell', ...
            'Unable to parse cell input to set BC.');
    end
    
else
    error('CHEBFUN:chebop:parsebc:unknown', 'Unsupported format of BCs.')

end

end
