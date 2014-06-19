classdef (InferiorClasses = {?chebfun, ?operatorBlock, ?functionalBlock}) chebmatrix
%CHEBMATRIX   Compound matrix for operators, CHEBFUNs, and scalars.
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
%   documentation of CHEBMATRIX/CELLFUN().
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
% DEVELOPER NOTES:
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS PROPERTIES:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties ( Access = public )
        % A cell used to store the components of a CHEBMATRIX internally.
        blocks = {}   % ( Cell Array )
        % Domain of the CHEBMATRIX.
        domain = []   % ( Kx1 double )
    end
    
    properties ( Dependent = true )
        % DIFFORDER is a dependent property.
        diffOrder
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        % Constructor.
        function A = chebmatrix(data, dom)
            
            data = parseData(data);
            
            if ( isempty(data) )
                % Return an empty CHEBMATRIX:
                return
                
            elseif ( iscell(data) )
                % E.g., CHEBMATRIX({A, f ; g.' m});
                A.blocks = data;
                A.domain = A.mergeDomains(data{:});
                
            elseif ( isa(data, 'chebmatrix') )
                % Simply return the same CHEBMATRIX:
                A.blocks = data.blocks;
                A.domain = data.domain;
                % NOTE: We can't set A = data because this breaks subclass
                % overloads of the constructor.
                
            else
                % E.g., CHEBMATRIX(f), CHEBMATRIX(a), or CHEBMATRIX(L).
                A.blocks = {data};
                try 
                    A.domain = data.domain;
                catch
                    % Forgive the fact that data doesn't have a domain.
                end
                
            end

            if ( nargin > 1 )
                A.domain = A.mergeDomains(A, dom);
            end
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS IMPLEMENTED IN THIS FILE:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
                   
        function A = set.domain(A, d)
        %SET.DOMAIN   Insert breakpoints in the domain of the CHEBMATRIX.
        %   We don't allow removing breakpoints, or changing endpoints.
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
            % Loop through all elements.
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
                error('CHEBFUN:CHEBMATRIX:end:ngt2', ...
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
        %   tranposed (which is consistent with the built in CELL class).
            for k = 1:numel(A.blocks)
                % TODO: TRANSPOSE() or CTRANSPOSE()?
                A.blocks{k} = transpose(A.blocks{k});
            end
            A.blocks = transpose(A.blocks);
        end
        
        function varargout = loglog(A, varargin)
        %LOGLOG   Log-log plot of a CHEBMATRIX.
        %   Note that no warning is thrown for negative data.

            s = cellfun(@(b) min(size(b)), A.blocks);
            % Throw an error if A contains linear operators.
            if ( ~all(isfinite(s(:))) )
                error('CHEBFUN:CHEBMATRIX:loglog:infBlock', ...
                        'LOGLOG plot of infinite blocks is not supported.')
            end
            
            % Standard CHEBMATRIX/PLOT():
            [h{1:3}] = plot(A, varargin{:});
            % Strip negative data:
            for j = 1:3 
                for k = 1:numel(h{j})
                    xData = get(h{j}(k), 'xData');
                    xData(xData < 0) = 0;
                    set(h{j}(k), 'yData', xData);
                    yData = get(h{j}(k), 'yData');
                    yData(yData < 0) = 0;
                    set(h{j}(k), 'yData', yData);
                end
            end
            % Set the X and Y Scale to be logarithmic:
            set(gca, 'XScale', 'log', 'YScale', 'log');
            if ( nargout > 0 )
                varargout = h;
            end 

        end
        
        function out = num2cell(A)
        %NUM2CELL   Extract block entries of a CHEBMATRIX:
            out = A.blocks;
        end
        
        function varargout = semilogx(A, varargin)
        %SEMILOGX   Semilogx plot of a CHEBMATRIX.
        %   Note that no warning is thrown for negative data.
            
            s = cellfun(@(b) min(size(b)), A.blocks);
            % Throw an error if A contains linear operators.
            if ( ~all(isfinite(s(:))) )
                error('CHEBFUN:CHEBMATRIX:semilogx:infBlock', ...
                        'SEMILOGX plot of infinite blocks is not supported.')
            end

            % Standard CHEBMATRIX/PLOT():
            [h{1:3}] = plot(A, varargin{:});
            % Strip negative data:
            for j = 1:3 % TODO: This
                for k = 1:numel(h{j})
                    xData = get(h{j}(k), 'xData');
                    xData(xData < 0) = 0;
                    set(h{j}(k), 'xData', xData);
                end
            end
            % Set the XScale to be logarithmic:
            set(gca, 'XScale', 'log');
            if ( nargout > 0 )
                varargout = h;
            end
                
        end
        
        function varargout = semilogy(A, varargin)
        %SEMILOGY   Semilogy plot of a CHEBMATRIX.
        %   Note that no warning is thrown for negative data.
        
            s = cellfun(@(b) min(size(b)), A.blocks);
            % Throw an error if A contains linear operators.
            if ( ~all(isfinite(s(:))) )
                error('CHEBFUN:CHEBMATRIX:semilogy:infBlock', ...
                        'SEMILOGY plot of infinite blocks is not supported.')
            end

            % Standard CHEBMATRIX/PLOT():
            [h{1:3}] = plot(A, varargin{:});
            % Strip negative data:
            for j = 1:3 % TODO: This
                for k = 1:numel(h{j})
                    yData = get(h{j}(k), 'yData');
                    yData(yData < 0) = 0;
                    set(h{j}(k), 'yData', yData);
                end
            end
            % Set the YScale to be logarithmic:
            set(gca, 'YScale', 'log');
            if ( nargout > 0 )
                varargout = h;
            end
        end

        
        function A = simplify(A, varargin)
        %SIMPLIFY   Simplify CHEBFUN components in a CHEBMATRIX.
        %   SIMPLIFY(A) simplifies the CHEBFUN elements in a CHEBMATRIX, as
        %   defined by CHEBFUN/SIMPLIFY(). Other entries in A are not affected.
        %
        % See also CHEBFUN/SIMPLIFY().
            isCheb = cellfun(@(v) isa(v, 'chebfun'), A.blocks);
            A.blocks(isCheb) = cellfun(@(v) simplify(v, varargin{:}), ...
                A.blocks(isCheb), 'UniformOutput', false);
        end
        
        function A = transpose(A)
        %TRANSPOSE   Transpose a CHEBMATRIX.
        %   TRANSPOSE(A) transposes A.BLOCKS, but _not_ its entries.
            A.blocks = ctranspose(A.blocks);
        end
        
        function out = vscale(A)
        %VSCALE   Vertical scale of components in a CHEBMATRIX.
        %   VSCALE(A) returns a matrix containing the 'vertical scales' of the
        %   entries of a CHEBMATRIX as defined below:
        %       CHEBFUN: CHEBFUN/VSCALE()
        %       DOUBLE:  ABS()
        %       OPERATOR: NAN
            isDoub = cellfun(@(v) isfinite(max(size(v))), A.blocks);
            isCheb = cellfun(@(v) isfinite(min(size(v))), A.blocks) & ~isDoub;
            
            % TODO: Throw an error instead?
%             if ( ~all(isCheb | isDoub) )
%                 error('CHEBFUN:CHEBMATRIX:vscale:op', ...
%                     ['VSCALE is not defined for CHEBMATRIX objects ', ...
%                      'containing operators.']);
%             end
            
            out = NaN(size(A.blocks));
            out(isDoub) = cellfun(@(v) abs(v), A.blocks(isDoub));
            out(isCheb) = cellfun(@(v) vscale(v), A.blocks(isCheb));
        
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CELLFUN METHODS IMPLEMENTED IN THIS FILE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods 
        
        function A = cellfun(A, op)
        %CELLFUN   Apply an operation to each block of a CHEBMATRIX.
        %   CELLFUN(A, OP) applies the operator OP to each of the blocks of
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
            A = cellfun(A, @abs);
        end
        
        function A = acos(A)
            A = cellfun(A, @acos);
        end
        
        function A = acosd(A)
            A = cellfun(A, @acosd);
        end
        
        function A = acosh(A)
            A = cellfun(A, @acosh);
        end
        
        function A = acot(A)
            A = cellfun(A, @acot);
        end
        
        function A = acotd(A)
            A = cellfun(A, @acotd);
        end
        
        function A = acoth(A)
            A = cellfun(A, @acoth);
        end 
        
        function A = acsc(A)
            A = cellfun(A, @acsc);
        end
        
        function A = acscd(A)
            A = cellfun(A, @acscd);
        end
        
        function A = acsch(A)
            A = cellfun(A, @acsch);
        end
        
        function A = asec(A)
            A = cellfun(A, @asec);
        end
        
        function A = asecd(A)
            A = cellfun(A, @asecd);
        end
        
        function A = asech(A)
            A = cellfun(A, @asech);
        end
        
        function A = asin(A)
            A = cellfun(A, @asin);
        end
        
        function A = asind(A)
            A = cellfun(A, @asind);
        end
        
        function A = asinh(A)
            A = cellfun(A, @asinh);
        end
        
        function A = atan(A)
            A = cellfun(A, @atan);
        end
        
        function A = atand(A)
            A = cellfun(A, @atand);
        end
        
        function A = atanh(A)
            A = cellfun(A, @atanh);
        end
        
        function A = cos(A)
            A = cellfun(A, @cos);
        end 
        
        function A = cosd(A)
            A = cellfun(A, @cosd);
        end
        
        function A = cosh(A)
            A = cellfun(A, @cosh);
        end
        
        function A = cot(A)
            A = cellfun(A, @cot);
        end
        
        function A = cotd(A)
            A = cellfun(A, @cotd);
        end
        
        function A = coth(A)
            A = cellfun(A, @coth);
        end
        
        function A = csc(A)
            A = cellfun(A, @csc);
        end
        
        function A = cscd(A)
            A = cellfun(A, @cscd);
        end
        
        function A = csch(A)
            A = cellfun(A, @csch);
        end
        
        function A = diff(A, varargin)
            A = cellfun(A, @(A) diff(A, varargin{:}));
        end
        
        function A = erf(A)
            A = cellfun(A, @erf);
        end
        
        function A = erfc(A)
            A = cellfun(A, @erfc);
        end
        
        function A = erfcinv(A)
            A = cellfun(A, @erfcinv);
        end
        
        function A = erfcx(A)
            A = cellfun(A, @erfcx);
        end
        
        function A = erfinv(A)
            A = cellfun(A, @erfinv);
        end
        
        function A = exp(A)
            A = cellfun(A, @exp);
        end
        
        function A = expm1(A)
            A = cellfun(A, @expm1);
        end
        
        function A = fix(A)
            A = cellfun(A, @fix);
        end
        
        function A = floor(A)
            A = cellfun(A, @floor);
        end
        
        function A = heaviside(A)
            A = cellfun(A, @heaviside);
        end
        
        function A = imag(A)
            A = cellfun(A, @imag);
        end
        
        function A = log(A)
            A = cellfun(A, @log);
        end
        
        function A = log10(A)
            A = cellfun(A, @log10);
        end
        
        function A = log1p(A)
            A = cellfun(A, @log1p);
        end
        
        function A = log2(A)
            A = cellfun(A, @log2);
        end
        
        function A = max(A)
            A = cellfun(A, @max);
            if ( all(cellfun(@isnumeric, A.blocks)) )
                A = cell2mat(A.blocks);
            end
        end
        
        function A = min(A)
            A = cellfun(A, @min);
            if ( all(cellfun(@isnumeric, A.blocks)) )
                A = cell2mat(A.blocks);
            end
        end
        
        function A = mrdivide(A, b)
            if ( isnumeric(b) && isscalar(b) )
                A = cellfun(A, @(A) mrdivide(A, b));
            else
                error('CHEBFUN:CHEBMATRIX:mrdivide:notScalar', ...
                    'CHEBMATRIX/MRDIVIDE only supports division by scalars.');
            end
        end
        
        function A = power(A, b)
            if ( isnumeric(b) )
                A = cellfun(A, @(A) power(A, b));
            elseif ( isnumeric(A) )
                A = cellfun(b, @(b) power(A, b));
            else
                A.blocks = cellfun(@power, A.blocks, b.blocks, ...
                    'UniformOutput', false);
            end
        end
        
        function A = real(A)
            A = cellfun(A, @real);
        end
        
        function A = reallog(A)
            A = cellfun(A, @reallog);
        end
        
        function A = rdivide(A, b)
            if ( isnumeric(b) )
                if ( isscalar(b) )
                    A = cellfun(A, @(A) rdivide(A, b));
                else
                    if ( ~all(size(A.blocks) == size(b)) )
                        error('CHEBFUN:CHEBMATRIX:rdivide:dimensions', ...
                            'Matrix dimensions must agree.');
                    end
                    for k = 1:numel(b)
                        A.blocks{k} = A.blocks{k}/b(k);
                    end
                end    
            elseif ( isnumeric(A) )
                A = cellfun(b, @(b) rdivide(A, b));
            else
                A.blocks = cellfun(@rdivide, A.blocks, b.blocks, ...
                    'UniformOutput', false);
            end
        end
        
        function A = round(A)
            A = cellfun(A, @round);
        end
        
        function A = sec(A)
            A = cellfun(A, @sec);
        end
        
        function A = secd(A)
            A = cellfun(A, @secd);
        end
        
        function A = sech(A)
            A = cellfun(A, @sech);
        end
        
        function A = sign(A)
            A = cellfun(A, @sign);
        end
        
        function A = sin(A)
            A = cellfun(A, @sin);
        end
        
        function A = sinc(A)
            A = cellfun(A, @sinc);
        end
        
        function A = sind(A)
            A = cellfun(A, @sind);
        end
        
        function A = sinh(A)
            A = cellfun(A, @sinh);
        end
        
        function A = sqrt(A)
            A = cellfun(A, @sqrt);
        end
        
        function A = sum(A)
            A = cellfun(A, @sum);
        end
        
        function A = tan(A)
            A = cellfun(A, @tan);
        end
        
        function A = tand(A)
            A = cellfun(A, @tand);
        end
        
        function A = tanh(A)
            A = cellfun(A, @tanh);
        end
        
        function A = times(A, b)
            if ( isnumeric(b) )
                A = cellfun(A, @(A) times(A, b));
            elseif ( isnumeric(A) )
                A = cellfun(b, @(b) times(A, b));
            else
                A.blocks = cellfun(@times, A.blocks, b.blocks, ...
                    'UniformOutput', false);
            end
        end
        
        function A = uminus(A)
            A = cellfun(A, @uminus);
        end
        
        function A = uplus(A)
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS: 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = protected, Static = true )
        
        % Merges domains (union of breakpoints, while checking endpoints).
        d = mergeDomains(varargin)
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PRIVATE METHODS IMPLEMENTED IN THIS FILE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = parseData(data)
    if ( isa(data, 'chebmatrix') )
        % Do nothing.
    elseif ( ~iscell(data) )
        if ( isnumeric(data) || isa(data, 'chebfun') )
            data = num2cell(data);
        else
            % You're on your own, Chuck.
        end
    elseif ( ~all(cellfun(@numel, data) == 1) )
        error('CHEBFUN:CHEBMATRIX:input:nonscalarcell', ...
            'Entries in cell input to CHEBMATRIX constructor must be scalar.')
    else
        isCheb = cellfun(@(f) isa(f, 'chebfun'), data);
        if ( any( cellfun(@numColumns, data(isCheb)) > 1) )
            error('CHEBFUN:CHEBMATRIX:input:arraychebfun', ...
                ['CHEBFUN entries in cell input to CHEBMATRIX constructor ', ...
                 'must be scalar.'])
        end
    end
end
