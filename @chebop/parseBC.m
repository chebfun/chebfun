function result = parseBC(N, BC, type)
%PARSEBC  Parse boundary conditions for CHEBOP object.
%   This method is not intended for end users. For information about boundary
%   conditions in CHEBOPs, see CHEBOP.
%
%   This method is invoked by the set methods for LBC, RBC, and BC. The types
%   and meanings of the allowed input are described in the documentation for
%   CHEBOP. The result is either empty or a function handle that represents the
%   given condition as needed internally.

% Copyright 2017 by The University of Oxford and The Chebfun Developers. 
% See http://www.chebfun.org/ for Chebfun information.

numIn = nargin(N);

if ( isempty(BC) )
    result = [];
    
elseif ( isnumeric(BC) )
    % We reach here if we get assigned BCs like N.lbc = [1; 3].
    %   * In the scalar case, this means setting up conditions so that the
    %     solution equals the first entry of the vector at the endpoint, its
    %     derivative the second entry etc. In the example above, for a
    %     second order problem on [0,1], this amounts to
    %         u(0) = 1, u'(0) = 3.
    %   * In the system case, this means that the value of the first unknown 
    %     function equals the first entry of the vector, the second unknown
    %     function the second entry of the vector, etc. In the example above,
    %     for a coupled system on [0,1], this amounts to 
    %         u(0) = 1, v(0) = 3.
    %     Higher order coupled systems are not supported, as that would require
    %     inspecting the diffOrder of the linearized chebop, a (relatively)
    %     computionally expensive operation. That information only becomes
    %     available for free once we've called \ for starting to solve the
    %     problem.
    result = @(varargin) setupNumericalConditions(BC, varargin);
    
elseif ( isa(BC, 'function_handle') )
    % If we are dealing with a scalar problem where the independent variable is
    % not specified in the function handle arguments, allow also passing an
    % input function handle that takes one argument. Otherwise, we request that
    % the number of input to the BC function handle is one less than the number
    % of arguments to the OP part. A special case is when BCs are specified as a
    % vector, this leads to nargin(BC) = -1 (since we assigned the boundary
    % conditions as an anonymous function with a varargin argument), if this is
    % the case and the number of inputs for the BC don't match the operator, a
    % meaningful error will be thrown later on.
    if ( (numIn == 0) || ( (numIn == 1) && (nargin(BC) == 1) ) || ...
            ( strcmp(type,'lrbc') && (nargin(BC) == (numIn - 1)) ) || ...
            ( strcmp(type,'lrbc') && (nargin(BC) == -1) ) || ...
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

function out = setupNumericalConditions(BC, varargin)
% SETUPNUMERICALCONDITIONS  Return anonymous function for evaluating BCs.
%
% We need this function so that we can construct a function handle that allows
% specifying boundary conditions with a vector. In the scalar case, we will be
% subtracting the nth entry of the BC vector from the (n-1)st derivative of the
% solution. In the system case, we will be subtracting the nth entry of the
% vector from the nth unknown function at the endpoint:

u = varargin{1};
% Check if we're dealing with scalar or system case:
if (length(u) == 1)
    % Scalar case, subtract entries of BC from U and its derivatives:
    u = u{1};
    % We must always have at least one condition
    out = u - BC(1);
    for bcCounter = 2:length(BC)
        % Add to the boundary conditions vector. It's tricky to initialize the
        % OUT argument to the correct dimensions, as we need to be able to
        % evaluate this both with CHEBFUNs, ADCHEBFUNs and TREEVARs, hence a
        % cheeky growth of the array.
        out = [out; diff(u, bcCounter-1) - BC(bcCounter)]; %#ok<AGROW>
    end
else
    assert(length(u) == length(BC), ...
        'CHEBFUN:CHEBOP:parseBC:numberOfConditions', ...
        ['The number of initial/final conditions appears to be greater than' ...
        ' the number of variables in problem. Specifying boundary ' ...
        'conditions as a vector for systems is only supported for first ' ...
        'order systems.']);
    % We must always have at least one condition
    out = u{1} - BC(1);
    for bcCounter = 2:length(BC)
        % Add to the boundary conditions vector. Like above, it's tricky to
        % initialize the OUT argument to the correct dimensions, as we need to
        % be able to evaluate this both with CHEBFUNs, ADCHEBFUNs and TREEVARs.
        % So, allow ourselves to grow the array, by subtracting elements from
        % the BC vector from the solution components.
        out = [out; u{bcCounter} - BC(bcCounter)]; %#ok<AGROW>
    end
end

end
