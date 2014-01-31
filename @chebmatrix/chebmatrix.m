classdef (InferiorClasses = {?chebfun, ?operatorBlock, ?functionalBlock}) chebmatrix
%CHEBMATRIX Compound matrix for operators, chebfuns, and scalars.
%   A chebmatrix contains blocks that are linear operators, functionals,
%   chebfuns, or scalars. They are used to tie together multiple functions, or
%   functions and scalars, as unknowns in a system, and to express linear
%   operators on those objects.
%
%   Normally the CHEBMATRIX constructor is not called directly. Instead, one
%   uses the usual [ ] or concatenation commands familiar for matrices. Block
%   sizes must be compatible, where function dimensions have size Inf. All
%   functions and operators in a chebmatrix must share compatible domains; i.e.,
%   they should all have the same endpoints. The resulting chebmatrix domain
%   includes all the breakpoints of the constituent blocks.
% 
%   Chebmatrices obey the expected arithmetic operations, such as + and *, if
%   the sizes are appropriate. 
%
%   If a chebmatrix contains only chebfuns, then the CHEBFUN method converts it
%   to an array-valued chebfun so that plotting and other commands work on it
%   normally. 
%
%   Examples:
%     d = [-2 2];   % function domain
%     I = operatorBlock.eye(d);     % identity operator
%     D = operatorBlock.diff(d);    % differentiation operator
%     x = chebfun(@(x)x, d);        % the variable x on d
%     M = operatorBlock.mult(x.^2); % multiplication operator
%     S = functionalBlock.sum(d);   % integration functional 
%     E = functionalBlock.eval(d);  % evaluation functional generator
% 
%     u = [ exp(x); pi; sin(x) ];   %  function; scalar; function
%     A = [ I+D, abs(x), M;
%           S, 0, E(2);  
%           D, x.^2, I ];
%
%     dA = A.domain                 % includes breakpoint at zero
%     sz = size(A)                  % 3 by 3 block array
%     [Am, An] = blockSizes(A)      % sizes of the blocks
%
%     spy(A)                % show the block structures
%     matrix(A, [4 4])      % discretize with 4 points in each subdomain
%     matrix(A*u, [4 4]) 
%
%     A21 = A(2, 1);    % get just the (2,1) block
%     A21.domain        % no breakpoint
%     matrix(A21,6)     % Clenshaw-Curtis weights
%
%   See also CHEBOPPREF, LINOP, CHEBMATRIX.MATRIX, CHEBMATRIX.SPY.
    
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer notes
%
% CHEBMATRIX is the class that enables concatenating various classes of objects
% of Chebfun into a single object. It has four fields:
%   * blocks: A Matlab cell used to store the components.
%   * prefs: A CHEBOPPREF is attached to the object.
%   * domain: The domain of the componenents, including the union of all
%             breakpoints.
%   * diffOrder: A dependent property. diffOrder is matrix valued, with values
%                corresponding to the diffOrders of each components, so that
%                each order of differentiation gives +1, and each order of
%                anti-derivative gives -1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        % A cell used to store the components of a CHEBMATRIX internally.
        blocks = {}
        prefs = cheboppref 
        domain
    end
    
    properties (Dependent)
        % diffOrder is a dependent property
        diffOrder
    end
    
    methods
        
        % Constructor.
        function A = chebmatrix(data)
            if ( isempty(data) )
                return
            elseif ( iscell(data) )
                A.blocks = data;
                A.domain = A.mergeDomains(data{:});
            elseif ( isa(data, 'chebmatrix') )
                A.blocks = data.blocks;
                A.domain = data.domain;
            elseif ( isa(data, 'chebfun') || isa(data, 'linBlock') )
                A.blocks = {data};
                A.domain = data.domain;
            end
        end
            
        function A = set.domain(A, d)
            % We don't allow removing breakpoints, or changing endpoints.
            A.domain = A.mergeDomains(d, A.domain);
        end
        
        function d = get.diffOrder(L)
            d = getDiffOrder(L);
        end
               
        function d = getDiffOrder(A)
        % GETDIFFORDER      Differential order of each chebmatrix block. 
        %   Also accessible via property: A.diffOrder
            
            d = zeros(size(A));
            % Loop through all elements
            for j = 1:numel(A.blocks);
                if ( isa(A.blocks{j}, 'operatorBlock') )
                    d(j) = A.blocks{j}.diffOrder;
                end
            end
        end
        
        function out = isFunVariable(A, k)
        %ISFUNVARIABLE  Which variables of the chebmatrix are functions?
        %
        %   A chebmatrix can operate on other chebmatrices. Operator and
        %   function blocks are applied to function components, whereas
        %   functions and scalar blocks multiply scalar components. The output
        %   of this function is a logical vector that is 1 for those columns of
        %   the chebmatrix which expect to operate on function components, and 0
        %   for those that expect to multiply scalars.
            [rowSize, colSize] = blockSizes(A);
            out = isinf(colSize(1, :));
            if ( nargin > 1 )
                out = out(k);
            end
        end
        
    end
    
    methods (Static = true, Access = protected)
        
        % Merges domains (union of breakpoints, while checking endpoints)
        d = mergeDomains(varargin)
        
    end
end