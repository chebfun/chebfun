function pass = test_cheb2leg(pref)

if ( nargin == 0 ) 
    pref = chebfunpref();
end

% TODO: Test dimension?

% Choose a tolerance:
tol = 5e-12;

%% Small N:

% Test scalar conversion:
N = 20;
c_cheb = [zeros(N, 1) ; 1];
c_leg = chebtech.cheb2leg(c_cheb);
err = norm(c_cheb - c_leg, inf);
pass(1) = err < tol;

% Test an arbitrary vector against a stored value:
c_cheb = 1./(1:N).^2; 
c_cheb(2:2:end) = -c_cheb(2:2:end);
c_leg = chebtech.cheb2leg(c_cheb);
c_leg19 = -0.011460983274163;
err = abs(c_leg(19) - c_leg19)/abs(c_leg19);
pass(2) = err < tol;

% Test conversion back to cheb coeffs:
c_cheb2 = chebtech.leg2cheb(c_leg);
err = norm(c_cheb2.' - c_cheb, inf);
pass(3) = err < tol;

%% Large N:

N = 1000;
c_cheb = [zeros(N, 1) ; 1];
c_leg = chebtech.cheb2leg(c_cheb);
err = norm(c_cheb-c_leg, inf);
pass(4) = err < 10*tol;

% Test an arbitrary vector against a stored value:
c_cheb = 1./(1:N).^2; 
c_cheb(2:2:end) = -c_cheb(2:2:end);
c_leg = chebtech.cheb2leg(c_cheb);
c_leg442 = 8.239505429144573e-04;
err = abs(c_leg(442) - c_leg442)/abs(c_leg442);
pass(5) = err < tol;

% Test conversion back to cheb coeffs:
c_cheb2 = chebtech.leg2cheb(c_leg);
err = norm(c_cheb2.' - c_cheb, inf);
pass(6) = err < tol;

end
