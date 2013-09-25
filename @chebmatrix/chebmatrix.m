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
        
        %Signatures of externally defined methods.
        
        % Fundamental algebraic operations.
        C = mtimes(A,B)
        C = plus(A,B)
        C = uminus(A)      
        
        % Replace each block by its DIM-dimensional discretization.
        A = matrixBlocks(L,dim,dom,matrixType)
        
        % Construct a single matrix based on DIM-dimensional blocks. 
        A = matrix(L,dim,dom,matrixType)

        % Construct the discrete linear system at a particular dimension.
        [A,b,dom] = linSystem(L,f,dim,matrixType)
        
        % Solve a linear system (including chebfun constructions).
        u = linsolve(L,f,type)

        % Assign BCs at the left/right endpoints.
        L = lbc(L,f,value)
        L = rbc(L,f,value)
        
        % Concatenation
        C = horzcat(A,B)
        C = vertcat(A,B)
        %TODO: CAT()
        
        B = subsref(A,sr)
        
    end
    
    methods
        
        % Constructor.        
        function A = chebmatrix(blocks)
            A.fundomain = chebmatrix.mergeDomains(blocks);
            A.blocks = blocks;
        end
        
        function d = domain(L)
            % DOMAIN(L) returns the domain on which functions are defined for
            % the chebmatrix L.
            d = L.fundomain;
        end
        
        function k = numbc(L)
            % NUMBC(L) returns the number of constraints attached to L. 
            k = length(L.constraints);
        end
        
        function varargout = size(L,varargin)
            %SIZE Number of blocks within the chebmatrix.
            %
            % S = SIZE(L) returns both dimensions.
            % S = SIZE(L,K) returns Kth dimension (K=1,2). 
            % [M,N] = SIZE(L) returns both as scalars. 
           [varargout{1:nargout}] = size(L.blocks,varargin{:});
        end
        
        function varargout = blocksizes(A)
            %BLOCKSIZES Sizes of the blocks within the chebmatrix.
            %
            % BLOCKSIZES(L) returns a cell of 1x2 size vectors.
            % [M,N] = BLOCKSIZES(A) returns two matrices of row/column sizes.
            if nargout <= 1
                varargout = {cellfun(@size,A.blocks,'uniform',false)};
            else
                B = A.blocks;
                varargout{1} = cellfun(@(x)size(x,1),B(:,1));
                varargout{2} = cellfun(@(x)size(x,2),B(1,:));
            end
        end
               
        function display(A)
            [m,n] = size(A);
            fprintf('\n  %ix%i block chebmatrix of types:\n\n',m,n)
            disp( cellfun(@class,A.blocks,'uniform',false) )
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
         
        function C = minus(A,B)
            C = plus(A,-B);
        end
              
        function u = mldivide(L,f)
            u = linsolve(L,f);
        end
                      
        function L = bc(L,f,value)
            L.constraints(end+1) = struct('op',f,'value',value);
        end
                        
        function out = iszero(f)
            % TODO: Implement this properly (for linearity detection)
            out = false;
        end
        
    end
    
    methods (Access=private)
        
        % Multiply chebmatrix by scalar.
        C = scalartimes(A,z)
     
        % Find the differential orders of each equation (row).
        d = getEqnDiffOrders(L)
        
        % Find the differential orders of each variable (column).
        d = getVarDiffOrders(L)
        
        % Figure out how much to reduce dimension in each equation.
        d = getDownsampling(L)
        
        % Construct operators for generic continuity at each breakpoint.
        C = domainContinuity(L,maxorder)
 
        % Append proper breakpoint continuity conditions to a linear system. 
        L = appendContinuity(L)            

    end
    
    methods (Static,Access=private)
        
        % Union of all breakpoints, with "fuzzy" equality. 
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