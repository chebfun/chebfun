% Test file for SPINOP2:

function pass = test_spinop2()

% Construction from STRING for GL equation:
S = spinop2('GL');
N = S.nonlin;
pass(1) = strcmpi(func2str(N), '@(u)u-(1+1.5i)*u.*(abs(u).^2)');

% Construction from DOM/TSPAN:
dom = [0 2*pi 0 2*pi];
tspan = [0 1];
S = spinop2(dom, tspan);
pass(2) = isequal(S.domain, dom);
pass(3) = isequal(S.tspan, tspan);

% Test recursive SUBSREF:
pass(4) = isequal(S.domain([2 4]), [2*pi 2*pi]);

end
