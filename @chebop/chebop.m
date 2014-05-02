classdef (InferiorClasses = {?double}) chebop
%CHEBOP  CHEBOP class for representing operators on functions defined on [a,b].
% N = CHEBOP(OP) creates a CHEBOP object N with operator defined by OP, which
% should be a handle to a function (often created using an anonymous function)
% that accepts a chebfun or a chebmatrix consisting of chebfuns and scalars as
% input and returns a CHEBFUN (or CHEBMATRIX). The first input argument to OP is
% the independent variable X, while all others represent dependent functions of
% X; if only one input argument is accepted by OP, it is the dependent variable.
% In case of coupled systems, the function OP must return vertically, not
% horizontally, concatenated elements. Note, this differs from the V4 syntax.
%
% Examples:
%
%   @(x, u) x.*diff(u) + u;                             % one dependent variable
%   @(x, u, v, w) [ u.*diff(v) ; diff(u,2)+w; diff(w)-v ];    % 3 dependent vars
%   @(u) diff(u,2) - exp(u);                  % no explicit independent variable
%
% The number of elements in the output CHEBMATRIX should typically equal the
% number of dependent variables, whether specified as names or CHEBMATRIX
% elements (see section on parameter-dependent problems below). 
%
% By default, the operator acts on CHEBFUN objects defined on the domain [-1,1].
% CHEBOP(OP, D), for a vector D, gives a different domain, with breakpoints (if
% any) described by D.
%
% %% BOUNDARY CONDITIONS: %%
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
% elements (again, differing from V4 syntax).
%
% Examples:
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
%   @(x, u) [ u(1) - u(0) ; sum(x.*u) ] % set u(1) = u(0), and definite integral
%                                       % of xu over the whole interval = 0.
%
% BCs can also be assigned to the CHEBOP N after it has been constructed, by
% N.lbc = ..., N.rbc = ..., and N.bc = ... . This will overwrite the BCs
% currently stored in the field being assigned to, but not the other fields).
%
% CHEBOP(OP, ..., 'init', U) provides a CHEBFUN/CHEBMATRIX as a starting point
% for nonlinear iterations or a PDE solution. See CHEBOP/SOLVEBVP and
% CHEBOP/PDE15S for details.
%
% Note that many fields can be set after the chebop object N is created: N.op,
% N.lbc, N.rbc, N.bc, N.init can all be assigned new values. For example:
%
%    N = chebop(-5, 5);  % Constructs an empty chebop on the interval [-5,5]
%    N.op = @(x, u) 0.01*diff(u, 2) - x.*u;
%    N.bc = 'dirichlet';
%    plot(N\1)
%
% %% PARAMETER DEPENDENT PROBLEMS: %%
%
% TODO: Revisit help text on parameter problems.
%
% There is some support for solving systems of equations containing unknown
% parameters without the need to introduce extra equations into the system. 
%
% Example, y'' + x.*y + p = 0, y(-1) = 1, y'(-1) = 1, y(1) = 1 can be solved via
%    N = chebop(@(x, y, p) diff(y,2) + x.*y + p)
%    N.lbc = @(y, p) [y - 1 ; diff(y)];
%    N.rbc = @(y, p) y - 1;
%    plot(N\0)
%
% This syntax will work whenever p is not differentiated within N.op, i.e.,
% something like @(x,y,p) diff(p*diff(y)) will require a second equation
% explicitly enforcing that diff(p) = 0.
%
% See also CHEBOP/MTIMES, CHEBOP/MLDIVIDE, CHEBOPPREF.   

% Copyright 2014 by The University of Oxford and The Chebfun Developers. See
% http://www.chebfun.org/ for Chebfun information.

    properties ( GetAccess = 'public', SetAccess = 'public' )
        domain = [];    % Domain of the operator
        op = [];        % The operator
        lbc = [];       % Left boundary condition(s)
        rbc = [];       % Right boundary condition(s)
        bc = [];        % Other/internal/mixed boundary conditions
        init = [];      % Initial guess of a solution
        % Default discretization for linear(ized0 problems
        discretizationType = @colloc2;
    end
    
    %% CONSTRUCTOR:
    methods
        
        function N = chebop(op, dom, lbcIn, rbcIn, init)
        % CHEBOP constructor
        
            if ( nargin == 0 )
                return
            end
            
            % No domain passed:
            if ( nargin < 2 )
                if ( ~isnumeric(op) )
                    % Get default domain from CHEBPREF():
                    p = chebpref();
                    dom = p.domain;
                else
                    % DOM was passed, but no OP.
                    dom = op;
                    op = [];                    
                end
            elseif ( nargin == 2 && isnumeric(op) )
                dom = [op, dom];
                op = [];
            end
            
            % Assign operator and domain:
            N.op = op;
            N.domain = dom;
            
            % Assign BCs and INIT if they were passed:
            if ( nargin == 3 )
                % CHEBOP(OP, DOM, BC):
                N.bc = lbcIn;
            elseif ( nargin == 4 )
                if ( isa(rbcIn, 'function_handle') || ischar(rbcIn) || ...
                        isnumeric(rbcIn))
                    % CHEBOP(OP, DOM, LBC, RBC):
                    N.lbc = lbcIn;
                    N.rbc = rbcIn;
                else
                    % CHEBOP(OP, DOM, BC, INIT):
                    N.bc = lbcIn;
                    N.init = rbcIn;
                end
            elseif ( nargin >= 4 )
                % CHEBOP(OP, DOM, LBC, RBC, INIT):
                N.lbc = lbcIn;
                N.rbc = rbcIn;
                N.init = init;
            end
            
        end
    end
    
    %% METHODS IMPLEMENTED IN THIS FILE:
    
    methods
        
        function N = set.lbc(N, val)
        %CHEBOP.SET.LBC   Set left boundary condition of a CHEBOP.
        %   CHEBOP.SET.LBC offers more control of setting left boundary
        %   conditions than simply accessing the .lbc field, or using standard
        %   subsref.
        
            % Need to know the nargin of the CHEBOP:
            nIn = nargin(N);
            
            if ( isempty(val) )
                N.lbc = [];
                
            elseif ( isnumeric(val) )
                if ( nIn > 2 )
                    % Only allow scalar numerical values to be passed if we are
                    % dealing with a scalar problem.
                    error('CHEBFUN:CHEBOP:SETLBC', ...
                        'Can only assign scalar BCs to scalar problems');
                else
                    N.lbc = @(u) u - val;
                end
            
            elseif ( isa(val, 'function_handle') )
                % If we are dealing with a scalar problem where the independent
                % variable is not specified in the function handle arguments,
                % allow also passing an input function handle that takes one
                % argument. Otherwise, we request that the number of input to
                % the LBC function handle is one less than the number of
                % arguments to the OP part.
                if ( ( (nIn <= 1) && (nargin(val) == 1) ) || ...
                        (nargin(val) == (nIn - 1)) )
                    N.lbc = val;
                else
                    error('CHEBFUN:CHEBOP:SETLBC', ...
                    'Number of inputs to BCs do not match operator.');
                end
                
            elseif ( strcmpi(val, 'neumann') )
                if ( nIn <= 2 )
                    N.lbc = @(u) diff(u);
                else
                    error('CHEBFUN:CHEBOP:SETLBC', ...
                    'Can only assign scalar BCs to scalar problems.');
                end
                
            elseif ( strcmpi(val, 'dirichlet') )
                if ( nIn <= 2 )
                    N.lbc = @(u) u;
                else
                    error('CHEBFUN:CHEBOP:SETLBC', ...
                    'Can only assign scalar BCs to scalar problems.');
                end    
                
            else
                error('CHEBFUN:CHEBOP:SETLBC', 'Unsupported format of BCs.')
            end
            
        end
        
        function N = set.rbc(N, val)
        %CHEBOP.SET.RBC   Set rigt boundary condition of a CHEBOP.
        %   CHEBOP.SET.RBC offers more control of setting left boundary
        %   conditions than simply accessing the .rbc field, or using standard
        %   subsref.

            % Need to know the nargin of the CHEBOP:
            nIn = nargin(N);
            
            if ( isempty(val) )
                N.rbc = [];
                
            elseif ( isnumeric(val) )
                if ( nIn > 2 )
                    % Only allow scalar numerical values to be passed if we are
                    % dealing with a scalar problem.
                    error('CHEBFUN:CHEBOP:SETRBC', ...
                        'Can only assign scalar BCs to scalar problems');
                else
                    N.rbc = @(u) u - val;
                end
            
            elseif ( isa(val, 'function_handle') )
                % If we are dealing with a scalar problem where the independent
                % variable is not specified in the function handle arguments,
                % allow also passing an input function handle that takes one
                % argument. Otherwise, we request that the number of input to
                % the LBC function handle is one less than the number of
                % arguments to the OP part.
                if ( ( (nIn <= 1) && (nargin(val) == 1) ) || ...
                        (nargin(val) == (nIn - 1)) )
                    N.rbc = val;
                else
                    error('CHEBFUN:CHEBOP:SETRBC', ...
                    'Number of inputs to BCs do not match operator.');
                end
                
            elseif ( strcmpi(val, 'neumann') )
                if ( nIn <= 2 )
                    N.rbc = @(u) diff(u);
                else
                    error('CHEBFUN:CHEBOP:SETRBC', ...
                    'Can only assign scalar BCs to scalar problems.');
                end
                
            elseif ( strcmpi(val, 'dirichlet') )
                if ( nIn <= 2 )
                    N.rbc = @(u) u;
                else
                    error('CHEBFUN:CHEBOP:SETRBC', ...
                    'Can only assign scalar BCs to scalar problems.');
                end    
                
            else
                error('CHEBFUN:CHEBOP:SETRBC', 'Unsupported format of BCs.')
            end
        end
        
        function N = set.bc(N, val)
        %CHEBOP.SET.BC   Set right boundary condition of a CHEBOP.
        %   CHEBOP.SET.BC offers more control of setting right boundary
        %   conditions than simply accessing the .bc field, or using standard
        %   subsref.
        
            if ( isempty(val) )
                N.bc = [];
            
            elseif ( isnumeric(val) )
            % Allow passing numerical values to the .BC field, which will impose
            % the conditions at both the left and right end of the domain. This
            % is for backwards compatability.          
                N.lbc = val; %#ok<MCSUP>
                N.rbc = val; %#ok<MCSUP>
            
           elseif ( isa(val, 'function_handle') )
                if ( nargin(N) ~= nargin(val) )
                    % When assigning to the BC field, we request that the number
                    % of input arguments of the function handle is the same as
                    % the number of inputs to the OP field.
                    error('CHEBFUN:CHEBOP:SETBC', ...
                        'Number of inputs to BCs must match operator.');                    
                else
                    N.bc = val;
                end
            elseif ( strcmpi(val, 'periodic') )
                N.bc = 'periodic';    
            elseif ( strcmpi(val, 'dirichlet') || strcmpi(val, 'neumann') )
                N.bc = [];
                N.lbc = val; %#ok<MCSUP>
                N.rbc = val; %#ok<MCSUP>
            else
                error('CHEBFUN:CHEBOP:SETRBC', ...
                    'Unsupported format of BCs')
            end
        end   
        
        function N = set.op(N, val)
        %CHEBOP.SET.OP   Set the differential equation part of a CHEBOP.
        %   CHEBOP.SET.OP offers more control of setting the DE left boundary
        %   conditions than simply accessing the .op field, or using standard
        %   subsref.
        
            % We're happy with function handles
            if ( isa(val, 'function_handle') || isempty(val) )
                N.op = val;
            elseif ( iscell(val) )
                error('CHEBFUN:CHEBOP:setOp:type', ...
                    ['Specifying differential equation as a cell of ', ...
                     'anonymous \nfunctions is no longer supported. Please '...
                     'use the syntax \n   N.op = @(x,u,v) [diff(u,2) + v; ' ...
                     '...]\ninstead of \n   N.op = {@(x,u,v) diff(u) + v; ' ...
                     '@{x,u,v) ...}']);
            else
                error('CHEBFUN:CHEBOP:setOp:type', ...
                    'Unknown type of argument for .op field of a chebop.');
            end
        end   
        
        function N = set.init(N, val)
        %CHEBOP.SET.OP   Set the differential equation part of a CHEBOP.
        %   CHEBOP.SET.OP offers more control of setting the DE left boundary
        %   conditions than simply accessing the .op field, or using standard
        %   subsref.
        
            % We're happy with function handles
            if ( isa(val, 'chebfun') && size(val, 2) > 1 )
                val = chebmatrix(mat2cell(val).');
                warning('Please use vertical concatenation for initial guess.')
            end
            
            N.init = val;
                
        end   

        
    end
    
    %% METHODS IMPLEMENTED IN OTHER FILES:
    
    methods
        
        % Linearize a CHEBOP around a CHEBFUN u.
        [L, res, isLinear] = linearize(N, u, x, flag);  
        
        %\   Chebop backslash.
        u = mldivide(N, rhs, pref)
        
        % The number of input arguments to a CHEBOP .OP field.
        nIn = nargin(N)
        
    end
    
    %% STATIC HIDDEN METHODS:
        
    methods ( Static = true, Hidden = true )
        % TODO: These should be private methods as well
        
        % Convert RHS to a format used internally in chebop.
        newRHS = convertToRHS(rhs, residual)
        
        % Controls information displayed for Newton iterations
        [displayFig, displayTimer] = displayInfo(mode, varargin);
        
        % Display at the finish of Newton iteration.
        displayInfoFinal(u, delta, iterNo, errEstDE, errEstBC, displayFig, ...
            displayTimer, pref)
        
        % Display at the start of Newton iteration.        
        [displayFig, displayTimer] = displayInfoInit(u,pref);
        
        % Display during Newton iteration.        
        displayInfoIter(u, delta, iterNo, normdu, cFactor, errEst, lendu, ...
            lambda, lenu, displayFig, displayTimer, pref);
        
        % Display special information for linear problems
        displayInfoLinear(u, normRes, pref)
        
    end
        
    
    %% STATIC METHODS:
        
    methods ( Static = true ) 
        % TODO: These should be private methods as well
        
        % Solve a linear problem posed with CHEBOP.
        [u, info] = solvebvpLinear(L, rhs, residual, displayInfo, pref)
        
        
        function out = norm(f, type)   
        % Compute norm when using CHEBOP (useful because we don't have a norm
        % method for the CHEBMATRIX class)
        
        % F is probably a CHEBMATRIX (might in some cases be a CHEBFUN).
        % TYPE determines what norm we use (currently not in use).
        
        % TODO: This is probably no longer necessary?
        
            if ( isa(f, 'chebmatrix') )
                out = 0;
                for k = 1:numel(f.blocks)
                    if ( isnumeric(f.blocks{k}) )
                        out = max(out, f.blocks{k});
                    else
                        out = max(out, get(f.blocks{k}, 'vscale'));
                    end
                end
            else
                out = get(f, 'vscale');
            end
        end

    end
    
end

