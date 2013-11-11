classdef blockCoeff 
    properties ( Access = public )
        coeffs = [];
        domain = [-1, 1];
    end
    
    methods
        function A = blockCoeff(varargin)
            if isempty(varargin{1})
                % This call is used to create an empty object of the class,
                % so that its methods will be used to process the stack.
                return
            elseif isa(varargin{1},'linBlock')
                % Convert the given linBlock to its function form by
                % evaluating its stack.
                L = varargin{1};
                dummy = blockCoeff([]);
                dummy.domain = L.domain;
                A = L.stack( dummy );
            else
                f = varargin{1};
                if ( ~iscell(f) )
                    f = {f};
                end
                % Called with data. Create a regular object. 
                A.coeffs = f;
                A.domain = varargin{2};
            end
        end

 
        % TODO: Domains need to be merged, or read off from the coeffs.
        % TODO: Use quasimatrix or chebmatrix
        
        function C = mtimes(A, B)
            if ( isempty(A.coeffs) || isempty(B.coeffs) )
                C = blockCoeff([]);
                return
            end
            
            c = B.coeffs;
            for k = 1:numel(B.coeffs)
                c{k} = A.coeffs{end}.*c{k};
            end
            % CHECK: Can the coeffs property ever be empty?
            z = { chebfun(0, A.coeffs{1}.domain) };
            for j = 1:numel(A.coeffs)-1
                B = diff(B);
                c = [z, c];
                for k = 1:numel(B.coeffs)
                    c{k} = c{k} + A.coeffs{end-j}.*B.coeffs{k};
                end
            end
            
            C = blockCoeff(c, A.domain);
        end
        
        function C = plus(A, B)
            sA = numel(A.coeffs);
            sB = numel(B.coeffs);

            if ( (sA == 0) || (sB==0) )
                C = blockCoeff([]);
                return
            end
            
            z = { chebfun(0, A.coeffs{1}.domain) };
            if ( sA < sB )
                A.coeffs = [repmat(z, 1, sB-sA), A.coeffs];
                sA = sB;
            elseif ( sB < sA )
                B.coeffs = [repmat(z, 1, sA-sB), B.coeffs];
            end
            c = cell(1, sA);
            for k = 1:sA
                c{k} = A.coeffs{k} + B.coeffs{k};
            end
            C = blockCoeff(c, A.domain);
        end
        
        function A = uminus(A)
            for k = 1:numel(A.coeffs)
                A.coeffs{k} = -A.coeffs{k};
            end
        end

    end
    
    methods
        
        % These are the basic constructors. 
        
        function I = eye(A)
            I = blockCoeff( chebfun(1, A.domain), A.domain );
        end
        
        function I = zeros(A)
            I = blockCoeff( chebfun(0, A.domain), A.domain );
        end
        
        function F = mult(A,f)
            F = blockCoeff( f, A.domain );
        end
        
        function C = cumsum(A,m)
            % Not supported. 
            % TODO
            error('Not supported')
        end
        
        function D = diff(A,order)
            
            if ( nargin < 2 )
                order = 1;
            end
            
            if ( ~isempty(A.coeffs) )
                % This syntax is used by the mtimes method to apply
                % differentiation to an existing operator.
                D = A;
            else
                D = eye(A);
            end
                        
            % Differentiate the correct number of times.
            c = D.coeffs;
            for d = 1:order
                m = numel(c);
                c = c([1, 1:m]);
                for k = 2:m
                    c{k} = diff(c{k}) + c{k+1};
                end
                c{m+1} = diff(c{m+1});
            end
            D = blockCoeff(c, A.domain);
            
        end
        
        function S = sum(A)
            % Not supported. 
            error('Not supported')
        end
        
        function E = feval(A,location,direction)
            % Not supported. 
            error('Not supported')
        end
        
        function F = inner(f)
            % Not supported. 
            error('Not supported') 
        end


    end
    
end
        
