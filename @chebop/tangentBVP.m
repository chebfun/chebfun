function [t,tau] = tangentBVP(H, u, lambda, told, tauold, prefs)

% Find a tangent to the curve H(u,lambda)=0 at the given point, by solving
% a boundary-value problem.

d = domain(u);
x = chebfun(@(x) x, d);

S = linearize(H, [u; lambda], x);

% Create a chebmatrix for the constraint which goes at the bottom
Jm = [told'*diag(sum(d)) tauold];

% Store constraint
Scon = S.constraint;

% Derivative of augmented operator
S = linop([S;Jm]);

% Reassign constraint
Scon.values = 0*Scon.values;
S.constraint = Scon;

% Right hand side
rhs = [chebfun(0,d);1];

ttau = linsolve(S, rhs, prefs);

% Extract function and scalar
t = ttau{1};
tau = ttau{2};

scale = sqrt(t'*t + tau^2);

t = t/scale;
tau = tau/scale;

end