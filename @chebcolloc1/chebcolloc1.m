classdef chebcolloc1 < chebcolloc
%CHEBCOLLOC1    Collocation discretization on 1st kind points.
%   CHEBCOLLOC1 is an implementation of CHEBCOLLOC that implements spectral
%   collocation on 1st-kind Chebyshev points for differential and integral
%   operators.
%
%   Linear algebra operations generally take O(N^3) flops, where N is determined
%   automatically to resolve the solution. You can control the allowed values of
%   N through CHEBOPPREF.
%
% See also CHEBCOLLOC, COLLOC, OPDISCRETIZATION, CHEBOPPREF, CHEBOP.

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% No subclass-specific properties needed, and no special constructor either.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% CLASS CONSTRUCTOR:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = false )
        
        function disc = chebcolloc1(varargin)
            disc = disc@chebcolloc(varargin{:});
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% STATIC METHODS:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods ( Access = public, Static = true )
        
        function tech = returnTech()
            %RETURNTECH    Return the appropriate tech to use for CHEBCOLLOC1.
            tech = @chebtech1;
        end
        
        function D = diffmat(N, k)
            %DIFFMAT  Chebyshev differentiation matrix.
            %   D = DIFFMAT(N) is the matrix that maps function values at N
            %   Chebyshev points of the 1st kind to values of the derivative of
            %   the interpolating polynomial at those points.
            %
            %   D = DIFFMAT(N, K) is the same, but for the Kth derivative.
            %
            % See also COLLOC/BARYDIFFMAT.

            if ( nargin < 2 )
                k = 1;
            end

            x = chebtech1.chebpts(N);           % First kind points.
            w = chebtech1.barywts(N);           % Barycentric weights.
            t = chebtech1.angles(N);            % acos(x).
            D = chebcolloc.baryDiffMat(x, w, k, t); % Construct matrix.
            
        end
        
        function Q = cumsummat(N)
            %CUMSUMMAT  Chebyshev integration matrix.
            %   Q = CUMSUMMAT(N) is the matrix that maps function values at N
            %   Chebyshev points to values of the integral of the interpolating
            %   polynomial at those points, with the convention that the first
            %   value is zero.
            
            % [TODO]: More efficient implementation?
            % [TODO]: Implement this at the COLLOC level?
                        
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
            
            % Matrix mapping coeffs -> integral coeffs. Note that the highest
            % order term is truncated.
            k = 1:N;
            k2 = 2*(k-1);  k2(1) = 1;  % Avoid divide by zero
            B = diag(1./(2*k),-1) - diag(1./k2,1);
            v = ones(N,1); v(2:2:end) = -1;
            B(1,:) = sum( diag(v)*B(2:N+1,:), 1 );
            B(:,1) = 2*B(:,1);
            
            Q = T*B*Tinv;
            
            % Store:
            CACHE{N} = Q;
            
        end
        
    end
    
end
