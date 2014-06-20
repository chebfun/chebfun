classdef colloc1 < colloc
%COLLOC1    Collocation discretization on 1st kind points.
%   COLLOC1 is an implementation of COLLOC that implements spectral
%   collocation on 1st-kind Chebyshev points for differential and integral
%   operators.
%
%   Linear algebra operations generally take O(N^3) flops, where N is determined
%   automatically to resolve the solution. You can control the allowed values of
%   N through CHEBOPPREF.
%
% See also COLLOC, CHEBDISCRETIZATION, CHEBOPPREF, CHEBOP.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% No subclass-specific properties needed, and no special constructor either.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function disc = colloc1(varargin)
            disc = disc@colloc(varargin{:});
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        function D = diffmat(n, k)
            %DIFFMAT   Chebyshev differentiation matrix
            %   D = DIFFMAT(N) is the matrix that maps function values at N Chebyshev
            %   points to values of the derivative of the interpolating polynomial at
            %   those points.
            %
            %   D = DIFFMAT(N,K) is the same, but for the Kth derivative.

            % TODO: This can probably be factored out to the chebtech level.
            % TODO: Cache this?
            
            if ( nargin < 2 ), k = 1; end
            if ( n == 0 ), D = []; return, end
            if ( k == 0 )
                D = eye(n);
                return
            end
            if ( n == 1 ), D = 0; return, end
            
            % 1st-kind Chebyshev grid:
            x = chebtech1.chebpts(n);
            % 1st-kind Barycentric weights:
            w = chebtech1.barywts(n);
            
            ii = (1:n+1:n^2)';              % indices of diagonal
            Dx = bsxfun(@minus,x,x');       % all pairwise differences
            Dx(ii) = Dx(ii) + 1;            % add identity
            Dxi = 1./Dx;                    % reciprocal
            Dw = bsxfun(@rdivide,w.',w);    % pairwise divisions
            Dw(ii) = Dw(ii) - 1;            % subtract identity
            
            % k = 1
            D = Dw .* Dxi;
            D(ii) = 0; D(ii) = - sum(D,2);                  % negative sum trick
            
            if ( k > 1 )
                % 2nd order
                D = 2*D .* (repmat(D(ii),1,n) - Dxi);
                D(ii) = 0; D(ii) = - sum(D,2);              % negative sum trick
                
                % higher orders
                for m = 3:k
                    D = m*Dxi .* (Dw.*repmat(D(ii),1,n) - D);
                    D(ii) = 0; D(ii) = - sum(D,2);          % negative sum trick
                end
            end
            
        end
        
        function Q = cumsummat(N)
            %CUMSUMMAT  Chebyshev integration matrix.
            %   Q = CUMSUMMAT(N) is the matrix that maps function values at N Chebyshev
            %   points to values of the integral of the interpolating polynomial at those
            %   points, with the convention that the first value is zero.
            
            % TODO: More efficient implementation?
            % TODO: This is duplicated in a number of places.
            
            % Copyright 2014 by The University of Oxford and The Chebfun Developers.
            % See http://www.chebfun.org/ for Chebfun information.
            
            N = N-1;
            
            persistent CACHE  % Stores computed values for fast return
            if ( isempty(CACHE) ), CACHE = {}; end  % First call
            
            if ( length(CACHE) >= N && ~isempty(CACHE{N}) )
                Q = CACHE{N};
                return
            else
                CACHE{N} = [];
            end
            
            % Matrix mapping coeffs -> values.
            T = chebtech1.coeffs2vals(eye(N+1));
            
            % Matrix mapping values -> coeffs.
            Tinv = chebtech1.vals2coeffs(eye(N+1));
            
            % Matrix mapping coeffs -> integral coeffs. Note that the highest order term is
            % truncated.
            k = 1:N;
            k2 = 2*(k-1);  k2(1) = 1;  % Avoid divide by zero
            B = diag(1./(2*k),-1) - diag(1./k2,1);
            v = ones(N,1); v(2:2:end) = -1;
            B(1,:) = sum( diag(v)*B(2:N+1,:), 1 );
            B(:,1) = 2*B(:,1);
            
            Q = T*B(end:-1:1,end:-1:1)*Tinv;
            % Make exact:
            Q(1,:) = 0;
            
            %Store:
            CACHE{N} = Q;
            
        end
        
    end
    
end
