classdef (InferiorClasses = {?double}) chebop
%CHEBOP  CHEBOP class for representing operators on functions defined on [a,b].
%
% N = CHEBOP(OP) creates a chebop object N with operator defined by OP, which
% should be a handle to a function, (often created using an anonymous function),
% that accepts a chebfun or a chebmatrix consisting of chebfuns and scalars as
% input and returns a chebfun (or chebmatrix). The first input argument to OP is
% the independent variable X, while all others represent dependent functions of
% X; if only one input argument is accepted by OP, it is the dependent variable.
% In case of coupled systems, the function OP must return vertically, not
% horizontally, concatenated elements.
% Examples:
%
%   @(x, u) x.*diff(u) + u;                             % one dependent variable
%   @(x, u, v, w) [ u.*diff(v); diff(u,2)+w; diff(w)-v ];     % 3 dependent vars
%   @(x, u) [ u(:,1)-diff(u(:,2)); diff(u(:,1)).*u(:,2) ];          % chebmatrix
%   @(u) diff(u,2) - exp(u);                  % no explicit independent variable
%   @(u) [ diff(u(:,2))-u(:,1), diff(u(:,1)) ];        % no independent variable
%
% The number of elements in the output chebmatrix should equal the number of
% dependent variables, whether specified as names or chebmatrix elements. For
% systems of equations not given in chebmatrix form, the first input argument is
% always the independent variable.
%
% By default, the operator acts on chebfuns defined on the domain [-1,1].
% CHEBOP(OP, D), for a vector D, gives a different domain, with breakpoints (if
% any) described by D.
%
% CHEBOP(OP, D, LBC, RBC) specifies boundary conditions for functions at the
% left and right endpoints of the domain D. Possible values for LBC and RBC are:
%
%   []          : No condition (for only assigning LBC or RBC in constructor).
%   scalar      : All variables equal the given value.
%   'dirichlet' : All variables equal zero.
%   'neumann'   : All variables have derivative zero.
%   function    : A function handle that must accept all dependent variables as 
%                 given in OP and return a chebfun or chebmatrix. All elements
%                 of the result are evaluated at the endpoint, and for the
%                 solution of the BVP, they are made to equal zero.
%
% A boundary condition function may be nonlinear; it must not accept the
% independent variable X as an input. Again, in case of systems, the function
% describing the boundary conditions must return vertically concatenated
% elements. Some examples:
%
%   @(u) diff(u) - 2;               % set u' = 2 at the appropriate endpoint
%   @(u, v, w) [ u - 1; w ];        % set u = 1 and w = 0 at the endpoint
%   @(u) u(:,2) - diff(u(:,2));     % set u_2 - (u_2)' = 0 at the endpoint
%   @(u, v, w) u.*v - diff(w)       % set u*v=w' at the endpoint
%
% CHEBOP(OP, D, BC) gives boundary or other side conditions in an alternate
% form. Choices for BC are:
%
%   scalar      : All variables equal the given value at both endpoints.
%   'dirichlet' : All variables equal zero at both endpoints.
%   'neumann'   : All variables have derivative zero at both endpoints.
%   'periodic'  : Impose periodicity on all dependent variables.
%   function    : See below.
%
% Note that the 'dirichlet' and 'neumann' keywords impose behavior that may not
% be identical to the common understanding of Dirichlet or Neumann conditions in
% every problem.
%
% When BC is passed in the CHEBOP call, the more specialized fields LBC and RBC
% are ignored.
%
% If BC is given a function handle, then each condition must give points
% explicitly or otherwise evaluate to a scalar. The function handle must return
% a column vector, not a row vector. Example:
%    @(x, u) [ u(1) - u(0); sum(x.*u) ] % set u(1) = u(0), and definite integral
%                                         of xu over the whole interval = 0
%
% BCs can also be assigned to the chebop N after it has been constructed,
% by N.lbc = ..., N.rbc = ..., and N.bc = ... . This will overwrite the BCs
% currently stored in the field being assigned to, but not the other
% fields).
%
% CHEBOP(OP,...,'init',U) provides a chebfun/chebmatrix as a starting point for
% nonlinear iterations or a PDE solution. See CHEBOP/SOLVEBVP and CHEBOP/PDE15S
% for details.
%
% N = CHEBOP(..., 'dim', DIM) where DIM is an integer informs the chebop that it
% operates on chebmatrices of dimension "(Inf x Dim) x 1".
%
% Note that many fields can be set after the chebop object N is created:
% N.op, N.lbc, N.rbc, N.bc, N.init can all be assigned new values. For
% example:
%
%    N = chebop(-5,5);  % Constructs an empty chebop on the interval [-5,5]
%    N.op = @(x,u) 0.01*diff(u,2) - x.*u;
%    N.bc = 'dirichlet';
%    plot(N\1)
%
% TODO: Revisit help text on parameter problems.
%
% There is some support for solving systems of equations containing unknown
% parameters without the need to introduce extra equations into the system.
% For example, y''+x.*y+p = 0, y(-1) = 1, y'(-1) = 1, y(1) = 1 can be
% solved via
%    N = chebop(@(x,y,p)diff(y,2)+x.*y+p,@(y,p)[y-1,diff(y)],@(y,p)y-1);
%    plot(N\0)
% This syntax will work whenever p is not differentiated within N.op, i.e.,
% something like @(x,y,p) diff(p*diff(y)) will require a second equation
% explicitly enforcing that diff(p) = 0.
%
%
% See also chebop/mtimes, chebop/mldivide, cheboppref.   
    properties ( GetAccess = 'public', SetAccess = 'public' )
        domain = [];    % Domain of the operator
        op = [];        % The operator
        lbc = [];       % Left boundary condition(s)
        rbc = [];       % Right boundary condition(s)
        bc = [];        % Other/internal/mixed boundary conditions
        init = [];      % Initial guess of a solution
        % Default discretization for linear problems
        discretizationType = @colloc2;
    end
    
    methods
        
        function N = chebop(op, dom, bcIn, init)
        % CHEBOP constructor
            if ( nargin == 0 )
                return
            end
            
            % No domain passed.
            if ( nargin < 2 )
                % Need to access chebfunpref to create an operator on the
                % default domain if none is passed.
                p = chebpref();
                dom = p.domain;
            end
            
            % Assign operator and domain
            N.op = op;
            N.domain = dom;
            
            % Assign BCs if they were passed
            if ( nargin >= 3 )
                N.bc = bcIn;
            end
            
            % Assign initial guess of solution if it was passed
            if ( nargin >= 4 )
                N.init = init;
            end
            
        end
        
        function u = mldivide(N, rhs)
        %\      Chebop backslash.
        %
        % See also: chebop/solvebvp
            
            % Create a cheboppref...
            pref = cheboppref;
            % ... and pass it along the inputs to solvebvp
            u = solvebvp(N, rhs, pref);
        end
        
        function nin = nargin(N)
        % CHEBOP.nargin     The number of input arguments to a chebop
            nin = nargin(N.op);
        end
        
        function N = set.lbc(N, val)
        % CHEBOP.SET.LBC    
        %
        % Offers more control of setting left boundary conditions than simply
        % accessing the .lbc field, or using standard subsref.
        
            % Need to know the nargin of the CHEBOP
            nin = nargin(N);
            
            % Only allow scalar numerical values to be passed if we are dealing
            % with a scalar problem.
            if ( isnumeric(val) )
                if nin <= 2
                    N.lbc = @(u) u - val;
                else
                    error('CHEBFUN:CHEBOP:SETLBC', ...
                    'Can only assign scalar BCs to scalar problems');
                end
            
            % Otherwise the input better be a function handle.
            elseif ( isa(val,'function_handle') )
                % If we are dealing with a scalar problem where the independent
                % variable is not specified in the function handle arguments,
                % allow also passing an input function handle that takes one
                % argument. 
                %
                % Otherwise, we request that the number of input to the LBC
                % function handle is one less than the number of arguments to
                % the OP part.
                if ( ( nin == 1 && nargin(val) == 1) || ...
                        ( nargin(val) == nin - 1) )
                    N.lbc = val;
                else
                    error('CHEBFUN:CHEBOP:SETLBC', ...
                    'Number of inputs to BCs do not match operator.');
                end
            else
                error('CHEBFUN:CHEBOP:SETLBC', ...
                    'Unsupported format of BCs')
            end
        end
        
        function N = set.rbc(N, val)
        % CHEBOP.SET.RBC    
        %
        % Offers more control of setting right boundary conditions than simply
        % accessing the .rbc field, or using standard subsref.
        
            % Need to know the nargin of the CHEBOP
            nin = nargin(N);
            
            % Only allow scalar numerical values to be passed if we are dealing
            % with a scalar problem.
            if ( isnumeric(val) )
                if nin <= 2
                    N.rbc = @(u) u - val;
                else
                    error('CHEBFUN:CHEBOP:SETRBC', ...
                    'Can only assign scalar BCs to scalar problems');
                end
                
            % Otherwise the input better be a function handle.
            elseif ( isa(val,'function_handle') )
                % If we are dealing with a scalar problem where the independent
                % variable is not specified in the function handle arguments,
                % allow also passing an input function handle that takes one
                % argument. 
                %
                % Otherwise, we request that the number of input to the LBC
                % function handle is one less than the number of arguments to
                % the OP part.                
                if ( ( nin == 1 && nargin(val) == 1) || ...
                        ( nargin(val) == nin - 1) )
                    N.rbc = val;
                else
                    error('CHEBFUN:CHEBOP:SETRBC', ...
                    'Number of inputs to BCs do not match operator.');
                end
            else
                error('CHEBFUN:CHEBOP:SETRBC', ...
                    'Unsupported format of BCs')
            end
        end
        
        function N = set.bc(N, val)
        % CHEBOP.SET.BC    
        %
        % Offers more control of setting other conditions than simply accessing
        % the .bc field, or using standard subsref.
        
            % Don't allow passing numerical values to the .BC field.
            if ( isnumeric(val) )
                error('CHEBFUN:CHEBOP:SETBC', ...
                    'Can not assign numerical BCs to .bc field.');
            
            % When assigning to the BC field, we request that the number of
            % input arguments of the function handle is the same as the number
            % of inputs to the OP field.
            elseif ( isa(val, 'function_handle') )
                if ( nargin(N) == nargin(val) )
                    N.bc = val;
                else
                    error('CHEBFUN:CHEBOP:SETBC', ...
                    'Number of inputs to BCs must match operator.');
                end
            else
                error('CHEBFUN:CHEBOP:SETRBC', ...
                    'Unsupported format of BCs')
            end
        end
        
        % Linearise a CHEBOP around a CHEBFUN u.
        [L, res, isLinear] = linearise(N, x, u, flag);
        
        
    end
        
    methods (Static = true) % These should be private methods as well
        
        % Controls information displayed for Newton iterations
        [displayFig, displayTimer] = displayInfo(mode, varargin);
        
        % Display at the start of Newton iteration.        
        [displayFig, displayTimer] = displayInfoInit(u,pref);
        
        % Display during Newton iteration.        
        displayInfoIter(u, delta, iterNo, normdu, cFactor, errEst, lendu, ...
            lambda, lenu, displayFig, displayTimer, pref);

        
        % Display at the finish of Newton iteration.        
        displayInfoFinal(u, delta, iterNo, errEstDE, errEstBC, displayFig, ...
            displayTimer, pref)
        
        function newRHS = convertToRHS(rhs, residual)
            [numRow, numCol] = size(residual);
            
            
            if ( length(rhs) == 1 )
                % Allow a scalar RHS to be converted to a RHS of correct
                % dimensions
                rhs = repmat(rhs, numRow, numCol);
            elseif ~( all(size(rhs) == [numRow, numCol]) )
                error('CHEBFUN:CHEBOP:CONVERTTORHS', ...
                    'RHS does not match output dimensions of operator.');
            end
            
            rhsBlocks = cell(numRow, numCol);
            resBlocks = residual.blocks;
            
            dom = getDomain(residual);
            
            % Convert numerical values in RHS vector into chebmatrix
            for rhsCounter = 1:numRow
                if isa(resBlocks{rhsCounter}, 'chebfun')
                    % If corresponding block in the residual is a chebfun, the
                    % rhs must also be made to be a chebfun
                    rhsBlocks{rhsCounter} = chebfun(rhs(rhsCounter), dom);
                else
                    % Otherwise, the entry in the chebmatrix will be a scalar
                    rhsBlocks{rhsCounter} = rhs(rhsCounter);
                end
            end
            % Convert the cell-array to a chebmatrix
            newRHS = chebmatrix(rhsBlocks);
        end
    end
    
end

