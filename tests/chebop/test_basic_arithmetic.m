function pass = test_basic_arithmetic(pref)

% TODO: We currently only test scalar autonomous chebops.
%       We also do not check for errors boundary conditions are dealt with.

if ( nargin == 0 )
    pref = cheboppref();
end

dom = [0 3];
u = chebfun(@sin, dom);
A = chebop(@(u) u, dom);
B = chebop(@(u) u, dom);

% Basic plus:
C = A + B;
pass(1) = norm(A(u) + B(u) - C(u)) == 0;

% Basic minus:
C = A - B;
pass(2) = norm(C(u)) == 0;

% Basic uminus:
C = -B;
pass(3) = norm(B(u) + C(u)) == 0;

% Basic times (scalar):
C = 2*A;
pass(4) = norm(2*A(u) - C(u)) == 0;
C = 2.*A;
pass(5) = norm(2*A(u) - C(u)) == 0;

% Basic divide (scalar):
C = A/2;
pass(6) = norm(A(u)/2 - C(u)) == 0;
C = A./2;
pass(7) = norm(A(u)/2 - C(u)) == 0;

% eye
I = eye(A);
pass(8) = norm(u - I(u)) == 0;
pass(9) = norm((A-I)*u - (A(u)-u)) == 0;

end