function g = diff(f, varargin) 
%DIFF   Differentiation of a BALLFUN.
%   DIFF(F, DIM) computes the derivative of F. If DIM = 1, the
%   derivative is taken in the r-variable. If DIM = 2, the derivative
%   is taken in the lambda-variable. If DIM = 3, the derivative is
%   taken in the theta-variable.
%
%   F = DIFF( F, DIM, K) computes the kth derivatives of F in the variable
%   given by DIM.
%
%   F = DIFF( F, DIM, K, "cart") computes the kth derivatives of F in the 
%   cartesian variables x, y or z given by DIM.
%
% See also SUM, SUM2, SUM3.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Parse user inputs:
if ( nargin == 2 )
    dim = varargin{1};
    K = 1;
    Cart = false;
elseif (nargin == 3 && isnumeric(varargin{2}))
    dim = varargin{1};
    K = varargin{2};
    Cart = false;
elseif (nargin == 3 && ~isnumeric(varargin{2}))
    dim = varargin{1};
    K = 1;
    Cart = true;
else
    dim = varargin{1};
    K = varargin{2};
    Cart = true;
end

% Initialization:
F = f.coeffs;

% Implement higher derivatives as repeated (iterated) differentiation
for j = 1:K
    F = onediff(F, dim, Cart);
end

% Return the derivative ballfun function 
g = ballfun( F, 'coeffs' ); 
end

function F = onediff(F, dim, Cart)
    
    % Grab dimension of f: 
    [m, n, p] = size( F ); 

    if Cart == false
        % Matrix of derivative in the r direction 
        if dim == 1
            
            DC1 = ultraS.diffmat( m, 1); 
            S01 = ultraS.convertmat(m, 0, 0);
            % Compute:
            % for k = 1:p 
            %     F(:,:,k) = S01 \ ( DC1*F(:,:,k) ); 
            % end
            % In vectorized form this becomes: 
            F = reshape( S01 \ ( DC1*reshape(F, m, []) ), m, n, []);
            
        % Matrix of derivative in the lambda direction
        elseif dim == 2  
        %     DFd2lambda = (1i^d2*spdiags(((-floor(n/2):floor(n/2))').^d2, 0, n, n));
        %     for k = 1:p
        %        F(:,:,k) = F(:,:,k)*DFd2lambda; 
        %     end
            F = F.*repmat( 1i*((-floor(n/2):ceil(n/2)-1)), m, 1, p); 
        
        % Matrix of derivative in the theta direction 
        else
        %     DFd3theta = (1i^d3*spdiags(((-floor(p/2):floor(p/2))').^d3, 0, p, p)); 
        %     % Permute r and theta to do Matrix operations 
        %     F = permute(F,[3,2,1]); 
        %     for i = 1:m 
        %         % Derivative with respect to theta 
        %         F(:,:,i) = DFd3theta*F(:,:,i); 
        %     end 
        %     % Permute back 
        %     F = permute(F,[3,2,1]); 
            F = F.*repmat( reshape(1i.*(-floor(p/2):ceil(p/2)-1),[1 1 p]), m, n, 1);
        end
    else
        % If m and p are odd, make them even
        m_tilde = m + mod(m,2);
        n_tilde = n;
        p_tilde = p + mod(p,2);
        
        % The variable coefficients in the definitions of the derivatives means
        % that the length of the columns and rows will increase by one wave number
        % after taking the derivatives with respect to x, y and z. 
        m_tilde = m_tilde+2;
        n_tilde = n_tilde+2;
        p_tilde = p_tilde+2;
        
        % Expand F
        Fexp = zeros(m_tilde,n_tilde,p_tilde);
        Fexp(1:m,floor(n_tilde/2)+1-floor(n/2):floor(n_tilde/2)+n-floor(n/2), floor(p_tilde/2)+1-floor(p/2):floor(p_tilde/2)+p-floor(p/2)) = F;
        
        % Derivative wrt to r, lambda and theta
        dR = onediff(Fexp, 1, false);
        dLam = onediff(Fexp, 2, false);
        dTh = onediff(Fexp, 3, false);
        
        % Useful spectral matrices
        Mr = ultraS.multmat(m_tilde, [0;1], 0);
        MsinL = trigspec.multmat(n_tilde, [0.5i;0;-0.5i] ); 
        McosL = trigspec.multmat(n_tilde, [0.5;0;0.5] );
        MsinT = trigspec.multmat(p_tilde, [0.5i;0;-0.5i] ); 
        McosT = trigspec.multmat(p_tilde, [0.5;0;0.5] );
        
        % Matrix of derivative in the x direction 
        if dim == 1
            % Compute 
            % for k = 1:p_tilde 
            %    dR(:,:,k) = dR(:,:,k)*McosL.';
            %    dLam(:,:,k) = -Mr\dLam(:,:,k)*MsinL.';
            %    dTh(:,:,k) = Mr\dTh(:,:,k)*McosL.';
            % end
            % For speed, we do it like this: 
%             dR = permute(reshape(reshape(permute(dR,[1 3 2]),[],n_tilde)*McosL.',m_tilde,p_tilde,[]),[1 3 2]);
            dLam = -reshape( Mr \ reshape(dLam, m_tilde, []) , m_tilde, n_tilde, []);
            dTh = reshape( Mr \ reshape(dTh, m_tilde, []) , m_tilde, n_tilde, []);
            for k = 1:p_tilde 
                dR(:,:,k) = dR(:,:,k)*McosL;
                dLam(:,:,k) = dLam(:,:,k)*MsinL.';
                dTh(:,:,k) = dTh(:,:,k)*McosL.';
            end
            
            % Compute: 
            %dR = permute(dR,[3,2,1]);
            %dLam = permute(dLam,[3,2,1]);
            %dTh = permute(dTh,[3,2,1]);
            %
            % G = zeros(size(dR));
            % for k = 1:m_tilde
            %   G(:,:,k) = MsinT*dR(:,:,k) + MsinT\dLam(:,:,k) + McosT*dTh(:,:,k);
            % end
            % Permute back
            % Fexp = permute(G,[3,2,1]);
            % For speed, we compute it like this: 
            Fexp = reshape( reshape(dR, [], p_tilde)*MsinT.' +...
                            reshape(dLam, [], p_tilde)/MsinT.'+...
                            reshape(dTh, [], p_tilde)*McosT.', ...
                                                m_tilde, n_tilde, p_tilde);
            
        % Matrix of derivative in the y direction
        elseif dim == 2  
            
            % Compute: 
            % for k = 1:p_tilde
            %   dR(:,:,k) = dR(:,:,k)*MsinL.';
            %   dLam(:,:,k) = Mr\dLam(:,:,k)*McosL.';
            %   dTh(:,:,k) = Mr\dTh(:,:,k)*MsinL.';
            % end
            % For speed, we compute it like this: 
            for k = 1:p_tilde 
                dR(:,:,k) = dR(:,:,k)*MsinL.';
                dLam(:,:,k) = dLam(:,:,k)*McosL.';
                dTh(:,:,k) = dTh(:,:,k)*MsinL.';
            end
%             dR = permute(reshape(reshape(permute(dR,[1 3 2]),[],n_tilde)*MsinL.',m_tilde,p_tilde,[]),[1 3 2]);
            dLam = reshape( Mr \ reshape(dLam, m_tilde, []) , m_tilde, n_tilde, []);
%             dLam = permute(reshape(reshape(permute(dLam,[1 3 2]),[],n_tilde)*McosL.',m_tilde,p_tilde,[]),[1 3 2]);
            dTh = reshape( Mr \ reshape(dTh, m_tilde, []) , m_tilde, n_tilde, []);
%             dTh = permute(reshape(reshape(permute(dTh,[1 3 2]),[],n_tilde)*MsinL.',m_tilde,p_tilde,[]),[1 3 2]);        
            
            % Compute: 
            %  dR = permute(dR,[3,2,1]);
            %  dLam = permute(dLam,[3,2,1]);
            %  dTh = permute(dTh,[3,2,1]);
            %
            % G = zeros(size(dR));
            % for k = 1:m_tilde
            %   G(:,:,k) = MsinT*dR(:,:,k) + MsinT\dLam(:,:,k) + McosT*dTh(:,:,k);
            % end
            % Permute back
            % Fexp = permute(G,[3,2,1]);
            % For speed, we compute it like this: 
                             
            Fexp = reshape( reshape(dR, [], p_tilde)*MsinT.' +...
                            reshape(dLam, [], p_tilde)/MsinT.'+...
                            reshape(dTh, [], p_tilde)*McosT.', ...
                                                m_tilde, n_tilde, p_tilde);
            
        % Matrix of derivative in the z direction 
        else
            % Compute: 
            % dR = permute(dR,[1,3,2]);
            %  dTh = permute(dTh,[1,3,2]);
            %             
            % G = zeros(size(dR));
            % for k = 1:n_tilde
            % 	G(:,:,k) = dR(:,:,k)*McosT.' - Mr\dTh(:,:,k)*MsinT.';
            % end
            % Permute back
            % Fexp = permute(G,[1,3,2]);
            % For speed, we compute it like this: 
            dTh = reshape( Mr \ reshape(dTh, m_tilde, []) , m_tilde, n_tilde, []);
            
            Fexp = reshape( reshape(dR, [], p_tilde)*McosT.' -...
                   reshape(dTh, [], p_tilde)*MsinT.',...
                                             m_tilde, n_tilde, []);
        end
        
%         tol = eps;
%         
%         pmod = mod(p_tilde - p,2)+1; 
%         nmod = mod(n_tilde - n,2)+1;        
%         alias1 = Fexp(m+1:end,:,:); 
%         a1 = norm(alias1(:),'inf') > 10*tol;
%         alias2 = Fexp(:,2:end-nmod,:); 
%         a2 = norm(alias2(:),'inf') > 10*tol;
%         alias3 = Fexp(:,:,2:end-pmod); 
%         a3 = norm(alias3(:),'inf') > 10*tol;
% 
%         % Truncate Fexp
%         if ( a1 )
%             for k = 1:p_tilde
%                 Fexp(1:m,:,k) = chebtech2.alias(Fexp(:,:,k), m);
%             end
%         end
%         if ( a2 )
%             for k = 1:p_tilde 
%                 Fexp(:,2:end-nmod,k) = trigtech.alias(Fexp(:,:,k).',n).';
%             end
%         end
%         if ( a3 )
%             for j = 1:n
%                 vj = reshape(Fexp(:,j,:), m_tilde, p_tilde);
%                 vj = trigtech.alias(vj.', p).';
%                 Fexp(:,j,2:end-pmod) = reshape( vj, m_tilde, 1, p );
%             end
%         end
%         F = Fexp(1:m, 2:end-nmod, 2:end-pmod);
        F = Fexp;
    end
end