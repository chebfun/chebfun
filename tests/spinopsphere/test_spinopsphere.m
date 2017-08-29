% Test file for SPINOPSPHERE:

function pass = test_spinopsphere()

% Construction from STRING for GL equation:
S = spinopsphere('GL');
N = S.nonlin;
pass(1) = strcmpi(func2str(N), '@(u)u-(1+1.5i)*u.*(abs(u).^2)');

% Construction from TSPAN:
tspan = [0 1];
S = spinopsphere(tspan);
pass(2) = isequal(S.tspan, tspan);

end
