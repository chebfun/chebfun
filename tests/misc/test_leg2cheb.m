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
c_leg = [zeros(N, 1) ; 1]';
c_cheb = leg2cheb(c_leg);
err = norm(c_leg - c_cheb, inf);
pass(1) = err < tol;

% Test an arbitrary vector against a stored value:
c_leg = 1./(N:-1:1)'.^2; 
c_leg(2:2:end) = -c_leg(2:2:end);
c_cheb = leg2cheb(c_leg);
c_cheb2 = -0.087275909551917;
err = abs(c_cheb(2) - c_cheb2)/abs(c_cheb2);
pass(2) = err < tol;

% Test conversion back to cheb coeffs:
c_leg2 = cheb2leg(c_cheb);
err = norm(c_leg2 - c_leg, inf);
pass(3) = err < tol;

%% Large N:

N = 1000;
c_leg = [zeros(N, 1) ; 1]';
c_cheb = leg2cheb(c_leg);
err = norm(c_leg - c_cheb, inf);
pass(4) = err < 10*tol;

% Test an arbitrary vector against a stored value:
c_leg = 1./(N:-1:1)'.^2; 
c_leg(2:2:end) = -c_leg(2:2:end);
c_cheb = leg2cheb(c_leg);
c_cheb559 = 6.379508600687388e-04;
err = abs(c_cheb(559) - c_cheb559)/abs(c_cheb559);
pass(5) = err < tol;

% Test conversion back to cheb coeffs:
c_leg2 = cheb2leg(c_cheb);
err = norm(c_leg2 - c_leg, inf);
pass(6) = err < 10*tol;

% Test vectorization: 
A = rand(3,10); E = A; 
B = leg2cheb( A );
for jj = 1:size(A,2)
    E(:,jj) = leg2cheb( A(:, jj) ); 
end
pass(7) = ( norm( B - E ) < tol ); 

% Test vectorization: 
A = rand(514,10); E = A; 
B = leg2cheb( A );
for jj = 1:size(A,2)
    E(:,jj) = leg2cheb( A(:, jj) ); 
end
pass(8) = ( norm( B - E ) < tol ); 

%% Test normalization, small N: 
A = zeros(10, 1); A(1) = 1; 
B = leg2cheb( A, 1 ); 
f = chebfun( B , 'coeffs'); 
pass(9) = ( norm( f ) - 1 ) < tol; 

%% Test normalization, large N: 
A = zeros(1000, 1); A(1) = 1; 
B = leg2cheb( A, 'normalize' ); 
f = chebfun( B , 'coeffs'); 
pass(10) = ( norm( f ) - 1 ) < tol; 


%% Test normalization and inverse, small N: 
A = zeros(10,1); A(1) = 1; 
B = leg2cheb( A, 1); 
C = cheb2leg( B, 1); 
pass(11) = ( norm( A - C ) ) < tol; 

%% Test normalization and inverse, large N: 
A = zeros(1000,1); A(1) = 1; 
B = leg2cheb( A, 1); 
C = cheb2leg( B, 1); 
pass(11) = ( norm( A - C ) ) < 10*tol; 
end
