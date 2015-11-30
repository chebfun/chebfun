function result = parseBC(N, BC, type)
%PARSEBC  Parse boundary conditions for CHEBOP object.
%   This method is not intended for end users. For information about boundary
%   conditions in CHEBOPs, see CHEBOP.
%
%   This method is invoked by the set methods for LBC, RBC, and BC. The types
%   and meanings of the allowed input are described in the documentation for
%   CHEBOP. The result is either empty or a function handle that represents the
%   given condition as needed internally.

% Copyright 2015 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

numIn = nargin(N);

if ( isempty(BC) )
    result = [];
    
elseif ( isnumeric(BC) )
    % This means setting up condition so that the solution equals the first
    % entry of the vector at the endpoint, its derivative the second entry etc.
    if ( numIn > 2 )
        % Allow only if we are dealing with a scalar problem.
        error('CHEBFUN:CHEBOP:parseBC:numeric', ...
            'Can only assign scalar BCs to scalar problems');
    else
        result = @(u) setupNumericalConditions(u, BC);
    end
    
elseif ( isa(BC, 'function_handle') )
    % If we are dealing with a scalar problem where the independent variable is
    % not specified in the function handle arguments, allow also passing an
    % input function handle that takes one argument. Otherwise, we request that
    % the number of input to the BC function handle is one less than the number
    % of arguments to the OP part.
    if ( (numIn == 0) || ( (numIn == 1) && (nargin(BC) == 1) ) || ...
            ( strcmp(type,'lrbc') && (nargin(BC) == (numIn - 1)) ) || ...
            ( strcmp(type,'bc') && (nargin(BC) == numIn) ) || ...
            ( strcmp(type,'bc') && (nargin(BC) == numIn + 1) ) )
        % We've got an anonymous function we're happy with. Do we need to
        % vectorize it?
        if ( N.vectorize && ~strcmp(func2str(BC), vectorize(BC)) )
            result = N.vectorizeOp(BC);
        else
            result = BC;
        end
    else
        error('CHEBFUN:CHEBOP:parseBC:inputs', ...
            'Number of inputs to BCs do not match operator.');
    end
    
elseif ( strcmpi(BC, 'neumann') )
    % Homogeneous Neumann.
    if ( numIn <= 2 )
        result = @(u) diff(u);
    else
        error('CHEBFUN:CHEBOP:parseBC:neuman', ...
            'Can only assign scalar BCs to scalar problems.');
    end
    
elseif ( strcmpi(BC, 'dirichlet') )
    % Homogeneous Dirichlet. 
    if ( numIn <= 2 )
        result = @(u) u;
    else
        error('CHEBFUN:CHEBOP:parseBC:dirichlet', ...
            'Can only assign scalar BCs to scalar problems.');
    end
    
elseif ( iscell(BC) ) && ( length(BC) == 2 ) && ( ischar(BC{2}) ) && ...
        ( isnumeric(BC{1}) )
    % A cell may have a numerical value followed by a keyword of 'dirichlet' or
    % 'neumann'.
    
    % This behavior is retained only for backward compatibility and only for
    % problems with one variable. 
    warning('CHEBFUN:CHEBOP:parseBC:keywordbc',...
        ['Keyword/value specifications of boundary conditions are ', ...
         'deprecated and may be removed in future versions of Chebfun.'])
            
    if ( strcmpi(BC{2}, 'neumann') && ( numIn <= 2 ) )
        result = @(u) diff(u) - BC{1};
    elseif ( strcmpi(BC{2}, 'dirichlet') && ( numIn <= 2 ) )
        result = @(u) u - BC{1};
    else
         error('CHEBFUN:CHEBOP:parseBC:cell', ...
            'Unable to parse cell input to set BC.');
    end
    
else
    error('CHEBFUN:CHEBOP:parseBC:unknown', 'Unsupported format of BCs.')

end

end

function out = setupNumericalConditions(u, BC)
% SETUPNUMERICALCONDITIONS  Return anonymous function for evaluating BCs.
%
% We need this function so that we can construct a function handle that allows
% specifying boundary conditions with a vector. Essentially, we are subtracting
% the nth entry of the vector from the (n-1)st derivative of U.

% We must always have at least one condition
out = u - BC(1);
for bcCounter = 2:length(BC)
    % Add to the boundary conditions vector. It's tricky to initialize the OUT
    % argument to the correct dimensions, as we need to be able to evaluate this
    % both with CHEBFUNs, ADCHEBFUNs and TREEVARs, hence a cheeky growth of the
    % array.
    out = [out; diff(u, bcCounter-1) - BC(bcCounter)]; %#ok<AGROW>
end
end
