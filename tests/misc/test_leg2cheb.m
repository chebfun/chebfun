function pass = test_leg2cheb(pref)

if ( nargin == 0 ) 
    pref = chebfunpref();
end

% TODO: Test dimension?

% Chhose a tolerance:
tol = 1e-13;

%% Small N:

% Test scalar conversion:
N = 20;
c_leg = [zeros(N, 1) ; 1];
c_cheb = leg2cheb(c_leg);
err = norm(c_leg - c_cheb, inf);
pass(1) = err < tol;

% Test an arbitrary vector against a stored value:
c_leg = 1./(1:N).^2; 
c_leg(2:2:end) = -c_leg(2:2:end);
c_cheb = leg2cheb(c_leg);
c_cheb19 = 0.087275909551917;
err = abs(c_cheb(19) - c_cheb19)/abs(c_cheb19);
pass(2) = err < tol;

% Test conversion back to cheb coeffs:
c_leg2 = cheb2leg(c_cheb);
err = norm(c_leg2.' - c_leg, inf);
pass(3) = err < tol;

%% Large N:

N = 1000;
c_leg = [zeros(N, 1) ; 1];
c_cheb = leg2cheb(c_leg);
err = norm(c_leg - c_cheb, inf);
pass(4) = err < 10*tol;

% Test an arbitrary vector against a stored value:
c_leg = 1./(1:N).^2; 
c_leg(2:2:end) = -c_leg(2:2:end);
c_cheb = leg2cheb(c_leg);
c_cheb442 = -6.379508600687388e-04;
err = abs(c_cheb(442) - c_cheb442)/abs(c_cheb442);
pass(5) = err < tol;

% Test conversion back to cheb coeffs:
c_leg2 = cheb2leg(c_cheb);
err = norm(c_leg2.' - c_leg, inf);
pass(6) = err < 10*tol;

end
