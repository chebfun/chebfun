classdef (InferiorClasses = {?chebfun,?linopOperator,?linopFunctional}) chebmatrix
    % No size/compatability checking whatsoever!
    
    properties
        blocks = {};
        
        % The chebmatrix domain may contain breakpoints that no individual
        % block uses, because it merges the breakpoints.
        fundomain
        
        % Storage of boundary conditions or other conditions
        constraints = struct('op',{},'value',{});  % empty array of structs
        
    end
    
    methods
        function A = chebmatrix(blocks)
            A.fundomain = chebmatrix.mergeDomains(blocks);
            A.blocks = blocks;
        end
        
        function d = domain(A)
            d = A.fundomain;
        end
        
        function k = numbc(A)
            k = length(A.constraints);
        end
        
        function varargout = size(A,varargin)
            % S = SIZE(A)
            % S = SIZE(A,K)
            % [M,N] = SIZE(A)
            %
            % Returns the number of blocks in each or Kth dimension.
            [varargout{1:nargout}] = size(A.blocks,varargin{:});
        end
        
        function varargout = blocksizes(A)
            % S = BLOCKSIZES(A)       returns cell of 1x2 sizes
            % [M,N] = BLOCKSIZES(A)   returns vectors of row/column sizes
            if nargout <= 1
                varargout = {cellfun(@size,A.blocks,'uniform',false)};
            else
                B = A.blocks;
                varargout{1} = cellfun(@(x)size(x,1),B(:,1));
                varargout{2} = cellfun(@(x)size(x,2),B(1,:));
            end
        end
        
        function C = horzcat(A,B)
            if isa(A,'chebmatrix')
                b1 = A.blocks;
            else
                b1 = {A};
            end
            if isa(B,'chebmatrix')
                b2 = B.blocks;
            else
                b2 = {B};
            end
            C = chebmatrix( horzcat(b1,b2) );
        end
        
        function C = vertcat(A,B)
            C = chebmatrix( vertcat(A.blocks,B.blocks) );
        end
        
        function display(A)
            [m,n] = size(A);
            fprintf('\n  %ix%i block chebmatrix of types:\n\n',m,n)
            disp( cellfun(@class,A.blocks,'uniform',false) )
        end
        
        function B = subsref(A,sr)
            % A(I,J) returns the slice (submatrix) of A as with ordinary
            % matrices.
            switch(sr(1).type)
                case '()'
                    B = chebmatrix( subsref(A.blocks,sr(1)) );
                case '{}'
                    B = subsref(A.blocks,sr(1));
                otherwise
                    if strcmp(sr(1).subs,'blocks')
                        B = A.blocks;
                        if length(sr) > 1
                            B = subsref(B,sr(2));
                        end
                    end
            end
        end
        
        function C = mtimes(A,B)
            
            % Needs to be upgraded to understand block*chebmatrix as well.
            if isnumeric(A)  % suspect these aren't necessary...historical...
                C = scalartimes(B,A);
            elseif isnumeric(B)
                C = scalartimes(A,B);
            else
                if isa(A,'linop')
                    A = chebmatrix({A});
                elseif isa(B,'linop')
                    B = chebmatrix({B});
                end
                [m,n] = size(A);
                Adata = A.blocks;   % needed to avoid subsref call later
                Bdata = B.blocks;
                [n,p] = size(B);
                C = cell(m,p);
                for i = 1:m
                    for j = 1:p
                        % It's tricky to just start a sum with "zero",
                        % because we don't know if it's added to an operator or
                        % a functional.
                        u = Adata{i,1}*Bdata{1,j};
                        for k = 2:n
                            u = u + Adata{i,k}*Bdata{k,j};
                        end
                        C{i,j} = u;
                    end
                end
                C = chebmatrix(C);
            end
        end
        
        function C = plus(A,B)
            [m,n] = size(A);
            C = cell(m,n);
            for i = 1:m
                for j = 1:n
                    C{i,j} = A.blocks{i,j} + B.blocks{i,j};
                end
            end
            C = chebmatrix(C);
        end
        
        function C = uminus(A)
            [m,n] = size(A);
            C = cell(m,n);
            for i = 1:m
                for j = 1:n
                    C{i,j} = -A.blocks{i,j};
                end
            end
            C = chebmatrix(C);
        end
        
        function C = minus(A,B)
            C = plus(A,-B);
        end
        
        function L = matrixBlocks(A,dim,matrixType)
            % matrixBlocks(A,DIM) returns a cell array in which each block is represented by its
            % DIM-dimensional discretization.
            %
            % matrixBlocks(A,DIM,TYPE) uses the discretization type TYPE.
            
            data = A.blocks;
            if (nargin < 3)
                matrixType = linop.defaultDiscretization;
            end
            [m,n] = size(A);
            
            % Discretize each block in place.
            L = cell(m,n);
            for i = 1:m
                for j = 1:n
                    item = data{i,j};
                    if isa(item,'linop')
                        L{i,j} = matrix(item,dim,matrixType);
                    elseif isa(item,'chebfun')
                        x = chebpts(dim, item.domain);
                        L{i,j} = item(x);
                    else   % scalar
                        L{i,j} = item;
                    end
                end
            end
        end
        
        
        function A = matrix(L,dim,type)
            % MATRIX(A,DIM) returns a discretization of the chebmatrix A based on dimension DIM.
            %
            % MATRIX(A,DIM,TYPE) uses the discretization type TYPE.
            
            if nargin < 3
                type = linop.defaultDiscretization;
            end
            
            blocks = matrixBlocks(L,dim,type);
            
            try
                A = cell2mat(blocks);
            catch
                error('Discretization sizes are not compatible.')
            end
            
        end
        
        function [A,b] = linSystem(L,f,dim,type)
            if nargin < 4
                type = linop.defaultDiscretization;
            end
            
            if isa(f,'chebfun')
                f = chebmatrix({f});
            end
            
            L = appendContinuity(L);
            
            dom = domain(L);
            if (length(dim)==1) 
                dim = repmat(dim,1,length(dom)-1);
            end          
            
            Ablocks = matrixBlocks(L,dim,type);
            bblocks = matrixBlocks(f,dim,type);
            
            % Resize the blocks and append the constraints.
            [m,n] = size(L);
            rows = cell(m+numbc(L),1);
            
            % Resize the operator rows according to differential order.
            dummy = type([]);
            d = getDownsampling(L);
            for i = find( ~isnan(d) )
                M = cat(2,Ablocks{i,:},bblocks{i});
                rows{i} = dummy.resize( M, dim-d(i), dim, domain(L) );
            end
                        
            % Append the discrete constraints.
            bc = L.constraints;
            for i = 1:length(bc)
                M = [ matrix(bc(i).op,dim,type), bc(i).value ];
                rows{m+i} = M;
            end
            
            aug = cell2mat(rows);
            A = aug(:,1:end-1);
            b = aug(:,end);
            
        end
        
        function u = mldivide(L,f)
            u = solve(L,f);
        end
        
        function u = solve(L,f,type)
            if nargin < 3
                type = linop.defaultDiscretization;
            end
            
            dimVals = [32 64 128 256 360 512 720 1024];
            [rowSize,colSize] = blocksizes(L);
            isFunVariable = isinf(colSize);
            
            dom = domain(L);
            numint = length(dom)-1;
            
            for k = 1:length(dimVals)
                dim = dimVals(k);

                [A,b] = linSystem(L,f,dim,type);
                uDiscrete = A\b;
                
                % Break the discrete solution into chunks representing
                % functions and scalars. 
                m = colSize;
                m(isFunVariable) = dim*numint; % replace Inf with discrete size
                u = mat2cell(uDiscrete,m,1);

                uFun = u(isFunVariable);
                uVals = cell2mat(uFun');
                
                % Take an arbitrary linear combination.
                s = 1 ./ (3*(1:size(uVals,2)));
                v = uVals*s(:);
                
                obj = type([]);  % to get a static method
                isDone = true;
                epsLevel = 0;
                for i = 1:numint
                    vint = v( (i-1)*dim + (1:dim) );
                    [t1,t2] = obj.convergeTest(vint);
                    isDone = isDone && t1;
                    epsLevel = max(epsLevel,t2);
                end
                if isDone
                    break
                end
            end
            
            if ~isDone
                warning('Linear system solution may not have converged.')
            end
            
            % The variable u is a cell array with the different components of
            % the solution. Because each function component may be piecewise
            % defined, we will loop through one by one. 
            for k = find( isFunVariable )
                funvals = reshape( u{k}, dim, numint );
                f = chebfun( num2cell(funvals,1), dom );  % piecewise defined
                u{k} = simplify(f, epsLevel );
            end
             
            u = chebmatrix(u);
            
        end

        
        function output = spy(A)
            data = matrixBlocks(A,10);
            h = cellplot(data);
            
            % CELLPLOT seems to cover up the text representations of double
            % values. We give a positive z-value so that they sit on top again.
            % And we hide its big ugly box.
            for i = 2:length(h)
                if strcmp(get(h(i),'type'),'text')
                    set(h(i),'position',[0 0 1]+get(h(i),'position'))
                    set(h(i-1),'vis','off')
                end
            end
            
            if nargout > 0
                output = h;
            end
        end
        
        function L = bc(L,f,value)
            L.constraints(end+1) = struct('op',f,'value',value);
        end
        
        function L = lbc(L,f,value)
            d = domain(L);
            E = linop.evalAt(d(1),d);
            if nargin < 3
                value = f;
                f = E;
            else
                f = E*f;
            end
            L = L.bc(f,value);
        end
        
        function L = rbc(L,f,value)
            d = domain(L);
            E = linop.evalAt(d(end),d);
            if nargin < 3
                value = f;
                f = E;
            else
                f = E*f;
            end
            L = L.bc(f,value);
        end
        
    end
    
    methods (Access=private)
        
        function C = scalartimes(A,z)
            [m,n] = size(A);
            C = cell(m,n);
            for i = 1:m
                for j = 1:n
                    C{i,j} = z*A.blocks{i,j};
                end
            end
            C = chebmatrix(C);
        end
        
        function d = getRowDiffOrders(L)
            [m,n] = size(L);
            d = zeros(1,m);
            for i = 1:m
                for j = 1:n
                    block = L.blocks{i,j};
                    if isa(block,'linopOperator')
                        d(i) = max(d(i),block.diffOrder);
                    elseif isa(block,'linopFunctional') || isnumeric(block)
                        d(i) = NaN;
                    end
                end
            end
            
        end
        
        function d = getVarDiffOrders(L)
            [m,n] = size(L);
            d = zeros(1,n);
            for j = 1:n
                for i = 1:m
                    block = L.blocks{i,j};
                    if isa(block,'linopOperator')
                        d(j) = max(d(j),block.diffOrder);
                    elseif isa(block,'linopFunctional') || isnumeric(block)
                        d(i) = NaN;
                    end
                end
            end

        end
        
        function d = getDownsampling(L)
            [m,n] = size(L);
            dRow = getRowDiffOrders(L);
            dVar = getVarDiffOrders(L);
            % Each variable requires a downsampling contribution equal to its
            % differential order. The total diff. orders of the variables might
            % not equal the total of the rows. 
            totalRow = sum( dRow(~isnan(dRow)) );
            totalVar = sum( dVar(~isnan(dVar)) );
            if (totalRow == totalVar)
                d = dRow;
            else
                % The only reasonable thing is to spread out the D.O. reductions
                % as evenly as possible among the rows.
                d = NaN(1,m);
                isOp = ~isnan(dRow);
                numOp = sum(isOp);
                k = floor( totalVar / numOp );
                d(isOp) = k;  % even distribution of order
                i = rem(totalVar,numOp);  % leftover to be made up
                if i > 0
                    % Increment the first i rows
                    incr = find(isOp,i);
                    d(incr) = d(incr) + 1;
                end
            end
 
        end
        
        function C = continuity(L,maxorder)
            % Returns expressions of continuity conditions at
            % the breakpoints of the domain of L. 
            %   C{m,k} has the (m-1)th-order derivative at breakpoint k
            
            d = domain(L);
            A = linop.eye(d);
            D = linop.diff(d,1);
            for m = 0:maxorder
                for k = 2:length(d)-1
                    El = linop.evalAt(d(k),d,'-');
                    Er = linop.evalAt(d(k),d,'+');
                    C{m+1,k} = (El-Er)*A;
                end
                A = D*A;
            end
        end

        function L = appendContinuity(L)
            % Append smoothness constraints at domain breakpoints.
            d = getVarDiffOrders(L);
            dom = domain(L);
            if ( max(d) > 0 ) && ( length(dom) > 2 )
                C = continuity(L,max(d)-1);
                Z = linop.evalAt(dom(1),dom)*(linop.zeros(dom));
                zero = chebmatrix({Z});
                for var = 2:length(d)
                    zero = [zero,Z];
                end
                for var = 1:length(d)
                    B = zero;
                    for m = 0:d(var)-1
                        for k = 2:length(dom)-1
                            B.blocks{var} = C{m+1,k};
                            L = L.bc(B,0);
                        end
                    end
                end
            end
        end
            

    end
    
    methods (Static,Access=private)
        function d = mergeDomains(blocks)
            
            d = cellfun(@(A) getDomain(A),blocks,'uniform',false);
            
            function out = getDomain(A)
                if ( isnumeric(A) )
                    out = [NaN NaN];
                elseif ( isa(A, 'linop') )
                    out = A.fundomain;
                else
                    out = get(A, 'domain');
                end
            end
            
            % Collect the endpoints and take the outer hull.
            leftEnds = cellfun(@(x) x(1),d);
            left = min(leftEnds(:));
            rightEnds = cellfun(@(x) x(end),d);
            right = max(rightEnds(:));
            
            % We want to soften 'equality' relative to the domain length.
            tol = 100*eps*(right-left);
            
            % Extract all the interior breakpoints.
            d = cellfun(@(x) x(2:end-1),d,'uniform',false);
            
            % Find the unique ones (sorted).
            breakpoints = cat(2,d{:});
            breakpoints = unique(breakpoints);
            
            if ~isempty(breakpoints)
                % Remove all too close to the left endpoint.
                isClose = ( breakpoints - left < tol );
                breakpoints(isClose) = [];
                
                % Remove all too close to the right endpoint.
                isClose = ( right - breakpoints < tol );
                breakpoints(isClose) = [];
                
                % Remove interior points too close to one another.
                isClose =  diff(breakpoints < tol );
                breakpoints(isClose) = [];
            end
            
            % Put it all together.
            d = [left breakpoints right];
            
        end
    end
    
end