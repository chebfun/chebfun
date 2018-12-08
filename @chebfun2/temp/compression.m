
function [A, S, B] = compression(Z, d, Y, tol, rel)
% 
% X = Z*diag(d)*Y';  
%
% Compresses the factors of X:
% If rel = 1, a relative tolerance is used, 
% and ||X -A*S*B'||_2 \leq tol ||X||_2. 
%
% Otherwise,  ||X -A*S*B'||_2 \leq tol.

% written by Heather Wilber (heatherw3521@gmail.com)
% Jan. 2018

%share out diag if pivots are real

if any(~isreal(d))
    Z= Z*diag(d); 
else
ds = sign(d);
d = sqrt(abs(d));
Z = Z*diag(d); 
Y = (Y*diag(d.*ds));
end

%compression step
[QZ, RZ] = qr( Z, 0); 
[QY, RY] = qr( Y, 0);

% % XX * YY^T = (QX*RX) * (QY*RY)^T = QX*(RX*RY^T)*QY^T
inner_piece = RZ*RY'; 
[A, S, B ] = svd( inner_piece ); 
if rel ==1
    tol = tol*S(1,1); 
end
idx = find(diag(S)>tol, 1, 'last');
A = A(:,1:idx); 
B = B(:,1:idx); 
S = S(1:idx,1:idx); 
A = QZ * A;  
B = QY * B; 

end



