function pass = test_dst(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% TODO: include more extensive tests.

n = 5;
seedRNG(42)
R = rand(n);

tol = 50*eps;

pass(1) = norm(chebfun.idst(chebfun.dst(R,1),1) - R) < tol;
pass(2) = norm(chebfun.idst(chebfun.dst(R,2),2) - R) < tol;
pass(3) = norm(chebfun.idst(chebfun.dst(R,3),3) - R) < tol;
pass(4) = norm(chebfun.idst(chebfun.dst(R,4),4) - R) < tol;

N = 15; 
x = rand(N, 2); 
tol = 100*eps;

% DST-I
DST1 = sin(pi/(N+1)*((0:N-1)+1)'*((0:N-1)+1));
y = DST1 * x;
pass(5) = norm( y - chebfun.dst(x, 1) ) < tol; 

% DST-II
DST2 = sin(pi/N*((0:N-1)+1/2)' * ((0:N-1)+1));
y = DST2 * x; 
pass(6) = norm( y - chebfun.dst(x, 2) ) < tol; 

% DST-III
DST3 = sin(pi/N*((0:N-1)+1)'*((0:N-1)+1/2));
DST3(:, end) = .5*DST3(:, end); 
y = DST3 * x; 
pass(7) = norm( y - chebfun.dst(x, 3) ) < tol; 

% DST-IV
DST4 = sin(pi/N*((0:N-1)+1/2)'*((0:N-1)+1/2));
y = DST4 * x; 
% chebfun.dst(x, 4) 
% chebfun.dct(flipud(x), 4)

pass(8) = norm( y - chebfun.dst(x, 4) ) < 10*tol; 

% IDST-I 
y = DST1 \ x;
pass(9) = norm( y - chebfun.idst(x, 1) ) < tol; 

% Idst-II
y = DST2 \ x; 
pass(10) = norm( y - chebfun.idst(x, 2) ) < tol; 

% Idst-III
y = DST3 \ x; 
pass(11) = norm( y - chebfun.idst(x, 3) ) < tol; 

% Idst-IV
y = DST4 \ x; 
pass(12) = norm( y - chebfun.idst(x, 4) ) < tol; 

% Check complex inputs 
c = rand(10,1) + 1i*rand(10,1); 
for j = 1:4
    v1 = chebfun.dst(c, j);
    v2 = chebfun.dst(real(c),j) + 1i*chebfun.dst(imag(c),j);
    pass(12+j) = norm(v1 - v2) < tol; 
end
end
