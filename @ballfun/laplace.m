function g = laplace(f)
% LAPLACE Laplacian of a BALLFUN function
%   LAPLACE(f) is the laplacian of the BALLFUN function f

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

g = diff(f,1,2,"Cart")+diff(f,2,2,"Cart")+diff(f,3,2,"Cart");

%% Old code
% [m,n,p] = size(f);
% 
% % Expand m and p if they are odd
% m_tilde = m + mod(m,2);
% p_tilde = p + mod(p,2);
% 
% % Spectral matrices
% Mr = ultraS.multmat(m_tilde, [0;1], 0); 
% Mr2 = ultraS.multmat(m_tilde, [.5;0;.5], 0);
% DC1 = ultraS.diffmat( m_tilde, 1); 
% S01 = ultraS.convertmat(m_tilde, 0, 0);
% DC2 = ultraS.diffmat( m_tilde, 2); 
% S02 = ultraS.convertmat(m_tilde, 0, 1);
% DF1 = 1i*spdiags((-floor(p_tilde/2):floor(p_tilde/2))', 0, p_tilde, p_tilde);
% DF2 = -spdiags((-floor(p_tilde/2):floor(p_tilde/2)).^2', 0, p_tilde, p_tilde);
% Msin = trigspec.multmat(p_tilde, [.5i;0;-.5i]);
% Mcos = trigspec.multmat(p_tilde, [.5;0;.5]); 
% Msin2 = trigspec.multmat(p_tilde, [-.25;0;.5;0;-.25]);
% I = speye(p_tilde,p_tilde);
% 
% % Build the operator
% A = (S02\DC2) + 2*((S01*Mr)\DC1);
% 
% % Permute f to get f(r,th,lam)
% F = f.coeffs;
% F = permute(F,[1,3,2]);
% 
% % Expand f
% Fexp = zeros(m_tilde,p_tilde,n);
% Fexp(1:m,1+mod(p,2):end,:) = F;
% 
% % Create the output
% Gexp = zeros(m_tilde,p_tilde,n);
% 
% % Compute only the useful Fourier modes
% ListFourierLambda = [];
% for j = 1:n
%   if (max(max(abs(Fexp(:,:,j)))) > 1e-16)
%       ListFourierLambda = [ListFourierLambda j];
%   end
% end
% 
% % Loop over lambda
% for j = ListFourierLambda
%     % Fourier mode
%     k = j - floor(n/2) - 1;
%     Gexp(:,:,j) = A*Fexp(:,:,j) + Mr2\Fexp(:,:,j)*((DF1.'*Mcos.'/Msin.') + DF2.' - k^2*I/Msin2.');
% end
% 
% % Truncate
% G = Gexp(1:m,1+mod(p,2):end,:);
% 
% % Permute back
% G = permute(G,[1,3,2]);
% g = ballfun(G);
end
