% Test file for conformal2.m.

function pass = test_conformal2(pref)

if ( nargin < 1 )
    pref = chebfunpref();
end

% Save the current warning state and then clear it
[lastmsg, lastid] = lastwarn();
lastwarn('');

% First example from help text:
z = chebfun('exp(1i*pi*z)','trig');
C1 = z.*abs(1+.1*z^4); C2 = .5*z.*abs(1+.2*z^3);
[f,finv,rho] = conformal2(C1, C2);        

pass(1) = (abs(rho-.539197)<1e-3) && (abs(finv(f(1+0.1i))-(1+0.1i))<1e-3);

end
