function pass = test_iszero
% Checks whether iszero information for linear chebops is being passed and
% worked with properly.

% AB, 21/11/2010

d = domain(0, 1);

I = eye(d);
D = diff(d);
Z = zeros(d);

% Concatenations
A = [I D; I Z];
pass(1) = all(all(A.iszero == [0 0; 0 1]));

% Plus
B1 = I + Z;
pass(2) = ~B1.iszero;

B2 = [I Z] + [I D];
pass(3) = ~any(B2.iszero);

B3 = [Z Z] + [Z D];
pass(4) = all(B3.iszero == [1 0]);

% Times
C1 = I*[D Z];
pass(5) = all(C1.iszero == [0 1]);

C2 = [I;I]*Z;
pass(6) = all(C2.iszero == [1; 1]);

% Scalars
E1  = D*0;
pass(7) = E1.iszero;

E2 = Z + 2;
pass(8) = ~E2.iszero;

end
