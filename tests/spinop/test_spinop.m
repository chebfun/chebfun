% Test file for SPINOP/SPINOP:

function pass = test_spinop()

% Construction from STRING for KS equation:
S = spinop('ks');
L = S.linearPart;
pass(1) = strcmpi(func2str(L), '@(u)-diff(u,2)-diff(u,4)');

% Construction from DOM/TSPAN:
dom = [0 2*pi];
tspan = [0 1];
S = spinop(dom, tspan);
pass(2) = isequal(S.domain, dom);
pass(3) = isequal(S.tspan, tspan);

end
