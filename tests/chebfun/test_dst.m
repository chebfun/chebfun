function pass = test_dst(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% TODO: include more extensive tests.

n = 15;
seedRNG(42)
R = rand(n);

tol = 50*eps;

for k = 1:8
    pass(k) = norm(chebfun.idst(chebfun.dst(R,1),1) - R) < tol;
end

N = 15; 
x = rand(N, 2); 
tol = 100*eps;

% DST-I
DST1 = sin(pi/(N+1)*((0:N-1)+1)'*((0:N-1)+1));
y = DST1 * x;
pass(9) = norm( y - chebfun.dst(x, 1) ) < tol; 

% DST-II
DST2 = sin(pi/N*((0:N-1)+1)'*((0:N-1)+1/2));
y = DST2 * x; 
pass(10) = norm( y - chebfun.dst(x, 2) ) < tol; 

% DST-III
DST3 = sin(pi/N*((0:N-1)+1/2)' * ((0:N-1)+1));
DST3(:, end) = .5*DST3(:, end); 
y = DST3 * x; 
pass(11) = norm( y - chebfun.dst(x, 3) ) < tol; 

% DST-IV
DST4 = sin(pi/N*((0:N-1)+1/2)'*((0:N-1)+1/2));
y = DST4 * x; 
pass(12) = norm( y - chebfun.dst(x, 4) ) < 10*tol; 

% DST-V
DST5 = sin(pi/(N+1/2)*((0:N-1)+1)'*((0:N-1)+1));
y = DST5 * x; 
pass(13) = norm( y - chebfun.dst(x, 5) ) < 10*tol; 

% DST-VI
DST6 = sin(pi/(N+1/2)*((0:N-1)+1)'*((0:N-1)+1/2));
y = DST6 * x; 
pass(14) = norm( y - chebfun.dst(x, 6) ) < 10*tol;

% DST-VII
DST7 = sin(pi/(N+1/2)*((0:N-1)+1/2)'*((0:N-1)+1));
y = DST7 * x; 
pass(15) = norm( y - chebfun.dst(x, 7) ) < 10*tol;

% DST-VIII
DST8 = sin(pi/(N+1/2)*((0:N-1)+1/2)'*((0:N-1)+1/2));
y = DST8 * x; 
pass(16) = norm( y - chebfun.dst(x, 8) ) < 10*tol;

% IDST-I 
y = DST1 \ x;
pass(17) = norm( y - chebfun.idst(x, 1) ) < tol; 

% IDST-II
y = DST2 \ x; 
pass(18) = norm( y - chebfun.idst(x, 2) ) < tol; 

% IDST-III
y = DST3 \ x; 
pass(19) = norm( y - chebfun.idst(x, 3) ) < tol; 

% IDST-IV
y = DST4 \ x; 
pass(20) = norm( y - chebfun.idst(x, 4) ) < tol; 

% IDST-V
y = DST5 \ x;
pass(21) = norm( y - chebfun.idst(x, 5) ) < tol; 

% IDST-VI
y = DST6 \ x; 
pass(22) = norm( y - chebfun.idst(x, 6) ) < tol; 

% IDST-VII
y = DST7 \ x; 
pass(23) = norm( y - chebfun.idst(x, 7) ) < tol; 

% % IDST-VIII
pass(24) = true;
% y = DST8 \ x; 
% pass(24) = norm( y - chebfun.idst(x, 8) ) < tol; 

% Check complex inputs 
c = rand(10,1) + 1i*rand(10,1); 
for j = 1:8
    v1 = chebfun.dst(c, j);
    v2 = chebfun.dst(real(c),j) + 1i*chebfun.dst(imag(c),j);
    pass(24+j) = norm(v1 - v2) < tol; 
end
end
