classdef blockUS
    properties
        size = [];  % arbitrary, but fixed in any one instance
        domain = [-1 1];
    end
    
    methods
        function A = blockUS(varargin)
            % ultraspherical method for ChebT coefficients.
            if ( nargin > 1 )
                if isa(varargin{1},'linBlock')
                    L = varargin{1};
                    A.size = varargin{2};
                    A.domain = L.domain;
                    A = L.delayFun( A );  % not sure if this is needed...
                else
                    A.size = varargin{1};
                    A.domain = varargin{2};
                end
            end
        end
        
        function D = diff( A, K )
            %DIFFMAT(N,K,N), computes the kth order US derivative matrix
            
            n = sum(dim(A));
            if K > 0
                D = spdiags((0:n-1)',1,n,n);
                for s = 1:K-1
                    D = spdiags(2*s*ones(n,1),1,n,n)*D;
                end
            else
               D = speye(n,n); 
            end
            %             if nargin == 3
            %                 dom = varargin{1};
            %                 diffMat = ((2./diff(dom))^k)*diffMat;
            %             end
        end
        
        function S = convert( A, K1, K2 )
            %CONVERTMAT(A,K1,K2), convert C^(K1) to C^(K2)
            
            n = sum(dim(A));
            S = speye(n);
            for s = K1:K2
                S = spconvert(n,s) * S;
            end
        end
    end

    methods (Static)
        function B = resize(A, m, n)
            % chop off some rows and columns
            B = A(1:m,1:n);
        end
        function [isDone,epsLevel] = testConvergence(v)
            % TODO:
            
            isDone = 1; epsLevel = eps; 
        end
        function f_coeffs = discretizeFunction(f,dim,dom)
            if ( nargin < 3 )
                dom = f.domain;
            end
            f_coeffs = flipud(get(f,'coeffs'));
            % prolong/truncate. 
            if length(f_coeffs) > dim 
                f_coeffs = f_coeffs(1:dim); 
            else
                f_coeffs = [f_coeffs;zeros(dim-length(f_coeffs),1)];
            end
            
            % TODO: Mapping to correct US basis. 
        end
    end
end