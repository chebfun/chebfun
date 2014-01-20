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
%     I = linop.eye(d);   % identity operator block
%     D = linop.diff(d);  % differentiation operator block
%     x = chebfun(@(x)x,d);   % the variable x on d
%     M = linop.mult(x.^2);   % multiplication operator block
%     S = linop.sum(d);   % integration functional 
%     E = linop.eval(d);  % evaluation functional generator
% 
%     u = [ exp(x); pi; sin(x) ];   %  function; scalar; function
%     A = [ I+D, abs(x), M;
%           S, 0, E(2);  
%           D, x.^2, I ];
%
%     dA = A.domain             % includes breakpoint at zero
%     sz = size(A)              % 3 by 3 block array
%     [Am,An] = blockSizes(A)   % sizes of the blocks
%
%     spy(A)           % show the block structures
%     matrix(A,[4 4])  % discretize with 4 points in each subdomain
%     matrix(A*u,[4 4]) 
%
%     A21 = A(2,1);   % get just the (2,1) block
%     A21.domain      % no breakpoint
%     matrix(A21,6)   % Clenshaw-Curtis weights
%
%   See also CHEBOPPREF, LINOP, CHEBMATRIX.MATRIX, CHEBMATRIX.SPY.
    
%  Copyright 2013 by The University of Oxford and The Chebfun Developers.
%  See http://www.chebfun.org for Chebfun information.
    
    properties
        blocks = {};
        prefs = cheboppref;        
    end
    
    properties (Dependent)
        domain
        diffOrder
    end
    
    methods
        
        % Constructor.
        function A = chebmatrix(data)
            if isempty(data)
                return
            elseif isa(data,'chebmatrix')
                A.blocks = data.blocks;
            elseif isa(data,'chebfun') || isa(data,'linBlock')
                A.blocks = {data};
            elseif iscell(data)
                A.blocks = data;
            end
            
            % Run this to check domain compatability
            A.domain;
        end
        
        function d = get.domain(L)
            if ( isempty(L) )
                d = [];
            else
                isnum = cellfun(@isnumeric,L.blocks);
                blocks = L.blocks(~isnum);
                d = cellfun(@(x) x.domain,blocks,'uniform',false);
                d = chebfun.mergeDomains(d{:});
            end
        end
        
        function d = getDomain(L)
            % DOMAIN(L) returns the domain on which functions are defined for
            % the chebmatrix L.
            d = L.domain;
        end       
        
        function d = get.diffOrder(L)
            d = getDiffOrder(L);
        end
               
        function d = getDiffOrder(A)
            d = zeros(size(A));
            for j = 1:numel(A.blocks);
                if ( isa(A.blocks{j}, 'operatorBlock') )
                    d(j) = A.blocks{j}.diffOrder;
                end
            end
        end
        
        function out = isFunVariable(A,k)
            [rowSize, colSize] = blockSizes(A);
            out = isinf(colSize(1, :));
            if ( nargin > 1 )
                out = out(k);
            end
        end
        
    end
              
end