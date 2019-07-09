function X = adi( A, B, F, p, q )
%ADI   Alternating direction implicit method for solving AX-XB=F.
% X = ADI( A, B, F, p, q) solves the Sylvester equation
%
%            A*X - X*B = F
%
% using the ADI method with shift parameters p and q.

%% Reference: 
% [1] Lu, An, and Eugene L. Wachspress. 
%  "Solution of Lyapunov equations by alternating direction implicit iteration." 
%  Comp. & Math. with Appl., 21.9 (1991): pp. 43-58.

m = size(A, 1);
n = size(B, 1);
X = zeros(m, n);
Im = speye(m);
In = speye(n);
for j = 1:numel(p)
    X = (F-(A+q(j)*Im)*X) / (B+q(j)*In);
    X = (A+p(j)*Im) \ ( F - X*(B+p(j)*In) );
end

end
