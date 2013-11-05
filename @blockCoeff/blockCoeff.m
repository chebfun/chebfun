classdef blockCoeff 
    properties ( Access = public )
        coeffs = [];
        domain = [-1, 1];
    end
    
    methods
        function A = blockCoeff(f, domain)
            if ( isempty(f) )
                return
            elseif ( ~iscell(f) )
                f = {f};
            end
            A.coeffs = f;
            A.domain = domain;
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
                B = blockCoeff.diff(B);
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
    
    methods (Static)
        
        % These are the basic constructors. 
        
        function I = eye(domain)
            I = blockCoeff( chebfun(1, domain), domain );
        end
        
        function I = zeros(domain)
            I = blockCoeff( chebfun(0, domain), domain );
        end
        
        function F = mult(f)
            F = blockCoeff( f, f.domain );
        end
        
        function C = cumsum(domain,m)
            % Not supported. 
            C = blockCoeff([]);
        end
        
        function D = diff(domain,order)
            
            if ( nargin < 2 )
                order = 1;
            end
            
            if isa(domain,'blockCoeff')
                % This syntax is used by the mtimes method to apply
                % differentiation to an existing operator.
                D = domain;
            else
                D = blockCoeff.eye(domain);
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
            D = blockCoeff(c, domain);
            
        end
        
        function S = sum(domain)
            % Not supported. 
            S = blockCoeff([]);
        end
        
        function E = feval(domain,location,direction)
            % Not supported. 
            E = blockCoeff([]);
        end
        
        function F = inner(f)
            % Not supported. 
            F = blockCoeff([]); 
        end


    end
    
end
        
