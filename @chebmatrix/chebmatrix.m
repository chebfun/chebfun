classdef (InferiorClasses = {?chebfun, ?operatorBlock, ?functionalBlock}) chebmatrix
%CHEBMATRIX Compound matrix for operators, CHEBFUNs, and scalars.
%   A CHEBMATRIX contains blocks that are linear operators, functionals,
%   CHEBFUNs, or scalars. They are used to tie together multiple functions, or
%   functions and scalars, as unknowns in a system, and to express linear
%   operators on those objects.
%
%   Normally the CHEBMATRIX constructor is not called directly. Instead, one
%   uses the usual [ ] or concatenation commands familiar for matrices. Block
%   sizes must be compatible, where function dimensions have size Inf. All
%   functions and operators in a CHEBMATRIX must share compatible domains; i.e.,
%   they should all have the same endpoints. The resulting CHEBMATRIX domain
%   includes all the breakpoints of the constituent blocks.
% 
%   CHEBMATRIX object obey the expected arithmetic operations, such as + and *,
%   if the sizes are appropriate.
%
%   If a CHEBMATRIX contains only CHEBFUNs and doubles, then the CHEBFUN/PLOT
%   method converts it to an array-valued chebfun so that plotting. Some methods
%   can be applied directly to the blocks of A and return a CHEBMATRIX. See the
%   documentation of APPLY2BLOCKS().
%
%   Examples:
%     d = [-2, 2];                  % function domain
%     I = operatorBlock.eye(d);     % identity operator
%     D = operatorBlock.diff(d);    % differentiation operator
%     x = chebfun(@(x) x, d);       % the variable x on d
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
% See also CHEBOPPREF, LINOP, CHEBMATRIX.MATRIX, CHEBMATRIX.SPY.
    
% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Developer notes
%
% CHEBMATRIX is the class that enables concatenating various classes of objects
% of Chebfun into a single object. It has four fields:
%   * blocks:    A Matlab cell used to store the components.
%   * domain:    The domain of the componenents, including the union of all
%                breakpoints.
%   * diffOrder: A dependent property. diffOrder is matrix valued, with values
%                corresponding to the diffOrders of each components, so that
%                each order of differentiation gives +1, and each order of
%                anti-derivative gives -1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        % A cell used to store the components of a CHEBMATRIX internally.
        blocks = {}   % ( Cell Array )
        % Domain of the CHEBMATRIX
        domain = []   % ( Kx1 double )
    end
    
    properties ( Dependent )
        % DIFFORDER is a dependent property
        diffOrder
    end
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTOR % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        
        function A = chebmatrix(data)
            
            if ( isempty(data) )
                % Return an empty CHEBMATRIX:
                return
                
            elseif ( iscell(data) )
                % E.g., CHEBMATRIX({A, f ; g.' m});
                A.blocks = data;
                A.domain = A.mergeDomains(data{:});
                
            elseif ( isa(data, 'chebmatrix') )
                % Simply return the same matrix:
                A = data;
                
            else
                % E.g., CHEBMATRIX(f), CHEBMATRIX(a), or CHEBMATRIX(L).
                A.blocks = {data};
                try 
                    A.domain = data.domain;
                catch ME
                    % Forgive the fact that data doesn't have a domain.
                end
                
            end
            
        end
        
    end
    
%% %%%%%%%%%%%%%%%%%%%%% METHODS IMPLEMENTED IN THIS FILE %%%%%%%%%%%%%%%%%%%%%%

    methods
            
        function A = set.domain(A, d)
        %SET.DOMAIN   Insert breakpoints in the domain of the CHEBMATRIX:
            % We don't allow removing breakpoints, or changing endpoints.
            A.domain = A.mergeDomains(d, A.domain);
        end
        
        function d = get.diffOrder(L)
        %GET.DIFFORDER    Differential order of each CHEBMATRIX block. 
        %   Also accessible via property: get(A, 'diffOrder');
            d = getDiffOrder(L);
        end
               
        function d = getDiffOrder(A)
        %GETDIFFORDER    Differential order of each CHEBMATRIX block. 
        %   Also accessible via property: A.diffOrder;
            d = zeros(size(A));
            % Loop through all elements
            for j = 1:numel(A.blocks);
                if ( isa(A.blocks{j}, 'operatorBlock') )
                    d(j) = A.blocks{j}.diffOrder;
                end
            end
        end
        
        function e = end(A, k, n)
        %END   Last index of a CHEBMATRIX.
        %   END(A, K, N) is called for indexing expressions involving a
        %   CHEBMATRIX A when end is part of the K-th index out of N indices.
        %   For example, the expression A(end-1,:) calls A's end method with
        %   END(A, 1, 2). Note that N must be less than or equal to two.
            if ( n > 2 )
                error('CHEBFUN:end:ngt2', ...
                    'Index exceeds CHEBMATRIX dimensions.');
            end
            s = size(A.blocks);
            if ( n == 1 )
                e = prod(s);
            else
                e = s(k);
            end
        end
        
        function out = isFunVariable(A, k)
        %ISFUNVARIABLE   Which variables of the CHEBMATRIX are functions?
        %   A CHEBMATRIX can operate on other chebmatrices. Operator and
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
        
        function A = ctranspose(A)
        %CTRANSPOSE   Transpose a CHEBMATRIX.
        %   CTRANSPOSE(A) transposes A.BLOCKS and each of the entries in
        %   A.BLOCKS. Note the block entries are transposed, not conjugate
        %   tranposed.
            for k = 1:numel(A.blocks)
                % TODO: TRANSPOSE() or CTRANSPOSE()?
                A.blocks{k} = transpose(A.blocks{k});
            end
            A.blocks = transpose(A.blocks);
        end
        
        function A = transpose(A)
        %TRANSPOSE   Transpose a CHEBMATRIX.
        %   TRANSPOSE(A) transposes A.BLOCKS, but _not_ its entries.
            A.blocks = ctranspose(A.blocks);
        end
    end
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% BLOCK METHODS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods 
        
        function A = apply2blocks(A, op)
        %APPLY2BLOCKS   Apply an operation to each block of a CHEBMATRIX.
        %   APPLY2BLOCK(A, OP) applies the operator OP to each of the blocks of
        %   A. If OP is not defined for one of the block entry types, then an
        %   error is thrown.
        %
        %   The following methods are supported:
        %       ABS(), ACOS(), ACOSD(), ACOSH(), ACOT(), ACOTD(), ACOTH(),
        %       ACSC(), ACSD(), ACSH(), ASEC(), ASECD(), ASECH(), ASIN(),
        %       ASIND(), ASINH(), ATAN(), ATAND(), ATANH(), COS(), COSD(),
        %       COSH(), COT(), COTD(), COTH(), CSC(), CSCD(), CSCH(), DIFF(),
        %       ERF(), ERFC(), ERFCINV(), ERFCX(), ERFINV(), EXP(), EXPM1(),
        %       FIX(), FLOOR(), HEAVISIDE(), IMAG(), LOG(), LOG10(), LOG1P(),
        %       LOG2(), REAL(), REALLOG(), SEC(), SECD(), SECH(), SIGN(), SIN(),
        %       SINC(), SIND(), SINH(), SQRT(), SUM(), TAN(), TAND(), TANH(),
        %       UMINUS(), UPLUS().
        
            A.blocks = cellfun(op, A.blocks, 'UniformOutput', false);
            
        end
                
        function A = abs(A)
            A = apply2blocks(A, @abs);
        end
        function A = acos(A)
            A = apply2blocks(A, @acos);
        end
        function A = acosd(A)
            A = apply2blocks(A, @acosd);
        end
        function A = acosh(A)
            A = apply2blocks(A, @acosh);
        end
        function A = acot(A)
            A = apply2blocks(A, @acot);
        end
        function A = acotd(A)
            A = apply2blocks(A, @acotd);
        end
        function A = acoth(A)
            A = apply2blocks(A, @acoth);
        end 
        function A = acsc(A)
            A = apply2blocks(A, @acsc);
        end
        function A = acscd(A)
            A = apply2blocks(A, @acscd);
        end
        function A = acsch(A)
            A = apply2blocks(A, @acsch);
        end
        function A = asec(A)
            A = apply2blocks(A, @asec);
        end
        function A = asecd(A)
            A = apply2blocks(A, @asecd);
        end
        function A = asech(A)
            A = apply2blocks(A, @asech);
        end
        function A = asin(A)
            A = apply2blocks(A, @asin);
        end
        function A = asind(A)
            A = apply2blocks(A, @asind);
        end
        function A = asinh(A)
            A = apply2blocks(A, @asinh);
        end
        function A = atan(A)
            A = apply2blocks(A, @atan);
        end
        function A = atand(A)
            A = apply2blocks(A, @atand);
        end
        function A = atanh(A)
            A = apply2blocks(A, @atanh);
        end    
        function A = cos(A)
            A = apply2blocks(A, @cos);
        end    
        function A = cosd(A)
            A = apply2blocks(A, @cosd);
        end
        function A = cosh(A)
            A = apply2blocks(A, @cosh);
        end  
        function A = cot(A)
            A = apply2blocks(A, @cot);
        end    
        function A = cotd(A)
            A = apply2blocks(A, @cotd);
        end
        function A = coth(A)
            A = apply2blocks(A, @coth);
        end  
        function A = csc(A)
            A = apply2blocks(A, @csc);
        end    
        function A = cscd(A)
            A = apply2blocks(A, @cscd);
        end
        function A = csch(A)
            A = apply2blocks(A, @csch);
        end        
        function A = diff(A, varargin)
            A = apply2blocks(A, @(A) diff(A, varargin{:}));
        end  
        function A = erf(A)
            A = apply2blocks(A, @erf);
        end    
        function A = erfc(A)
            A = apply2blocks(A, @erfc);
        end
        function A = erfcinv(A)
            A = apply2blocks(A, @erfcinv);
        end 
        function A = erfcx(A)
            A = apply2blocks(A, @erfcx);
        end
        function A = erfinv(A)
            A = apply2blocks(A, @erfinv);
        end 
        function A = exp(A)
            A = apply2blocks(A, @exp);
        end
        function A = expm1(A)
            A = apply2blocks(A, @expm1);
        end 
        function A = fix(A)
            A = apply2blocks(A, @fix);
        end    
        function A = floor(A)
            A = apply2blocks(A, @floor);
        end         
        function A = heaviside(A)
            A = apply2blocks(A, @heaviside);
        end  
        function A = imag(A)
            A = apply2blocks(A, @imag);
        end
        function A = log(A)
            A = apply2blocks(A, @log);
        end         
        function A = log10(A)
            A = apply2blocks(A, @log10);
        end     
        function A = log1p(A)
            A = apply2blocks(A, @log1p);
        end     
        function A = log2(A)
            A = apply2blocks(A, @log2);
        end      
        function A = real(A)
            A = apply2blocks(A, @real);
        end
        function A = reallog(A)
            A = apply2blocks(A, @reallog);
        end        
        function A = round(A)
            A = apply2blocks(A, @round);
        end
        function A = sec(A)
            A = apply2blocks(A, @sec);
        end
        function A = secd(A)
            A = apply2blocks(A, @secd);
        end
        function A = sech(A)
            A = apply2blocks(A, @sech);
        end
        function A = sign(A)
            A = apply2blocks(A, @sign);
        end   
        function A = sin(A)
            A = apply2blocks(A, @sin);
        end
        function A = sinc(A)
            A = apply2blocks(A, @sinc);
        end
        function A = sind(A)
            A = apply2blocks(A, @sind);
        end
        function A = sinh(A)
            A = apply2blocks(A, @sinh);
        end
        function A = sqrt(A)
            A = apply2blocks(A, @sqrt);
        end        
        function A = sum(A)
            A = apply2blocks(A, @sum);
        end        
        function A = tan(A)
            A = apply2blocks(A, @tan);
        end
        function A = tand(A)
            A = apply2blocks(A, @tand);
        end
        function A = tanh(A)
            A = apply2blocks(A, @tanh);
        end
        function A = uminus(A)
            A = apply2blocks(A, @uminus);
        end       
        function A = uplus(A)
        end
        
    end
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% STATIC METHODS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    
    methods (Static = true, Access = protected)
        
        % Merges domains (union of breakpoints, while checking endpoints)
        d = mergeDomains(varargin)
        
    end
end
