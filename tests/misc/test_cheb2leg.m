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
c_cheb = [zeros(N, 1) ; 1]';
c_leg = cheb2leg(c_cheb);
err = norm(c_cheb - c_leg, inf);
pass(1) = err < tol;

% Test an arbitrary vector against a stored value:
c_cheb = 1./(N:-1:1)'.^2; 
c_cheb(2:2:end) = -c_cheb(2:2:end);
c_leg = cheb2leg(c_cheb);
c_leg2 = 0.011460983274163;
err = abs(c_leg(2) - c_leg2)/abs(c_leg2);
pass(2) = err < tol;

% Test conversion back to cheb coeffs:
c_cheb2 = leg2cheb(c_leg);
err = norm(c_cheb2 - c_cheb, inf);
pass(3) = err < tol;

%% Large N:

N = 1000;
c_cheb = [zeros(N, 1) ; 1]';
c_leg = cheb2leg(c_cheb);
err = norm(c_cheb-c_leg, inf);
pass(4) = err < 10*tol;

% Test an arbitrary vector against a stored value:
c_cheb = 1./(N:-1:1)'.^2; 
c_cheb(2:2:end) = -c_cheb(2:2:end);
c_leg = cheb2leg(c_cheb);
c_leg442 = -8.239505429144573e-04;
err = abs(c_leg(559) - c_leg442)/abs(c_leg442);
pass(5) = err < 10*tol;

% Test conversion back to cheb coeffs:
c_cheb2 = leg2cheb(c_leg);
err = norm(c_cheb2 - c_cheb, inf);
pass(6) = err < tol;

% Test vectorization: 
seedRNG(0)
A = rand(3,10); E = A; 
B = cheb2leg( A );
for jj = 1:size(A,2)
    E(:,jj) = cheb2leg( A(:, jj) ); 
end
pass(7) = ( norm( B - E ) < tol ); 

% Test vectorization: 
A = rand(514,10); E = A; 
B = cheb2leg( A );
for jj = 1:size(A,2)
    E(:,jj) = cheb2leg( A(:, jj) ); 
end
pass(8) = ( norm( B - E ) < tol ); 

% Test normalization: 
A = rand(10, 2); 
B = cheb2leg( A ); 
C = cheb2leg( A, 'norm'); 
D = cheb2leg( A, 'normalized'); 
pass(9) = ( norm( diag(1./(sqrt((0:9)' + 1/2)))*B - C ) < 10*tol ); 
pass(10) = ( norm( C - D ) < tol ); 

% Test normalization: 
A = rand(1000, 2);  
B = cheb2leg( A ); 
C = cheb2leg( A, 'norm'); 
D = cheb2leg( A, 'normalized'); 
pass(11) = ( norm( diag(1./(sqrt((0:999)' + 1/2)))*B - C ) < tol ); 
pass(12) = ( norm( C - D ) < tol ); 
end
