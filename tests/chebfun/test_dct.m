function pass = test_dct(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% TODO: include more extensive tests.

n = 5;
seedRNG(42)
R = rand(n);

tol = n*eps;

pass(1) = norm(chebfun.idct(chebfun.dct(R,1),1) - R) < tol;
pass(2) = norm(chebfun.idct(chebfun.dct(R,2),2) - R) < tol;
pass(3) = norm(chebfun.idct(chebfun.dct(R,3),3) - R) < tol;
pass(4) = norm(chebfun.idct(chebfun.dct(R,4),4) - R) < tol;

N = 13; 
x = rand(N, 2); 
tol = 50*eps;

% DCT-I
DCT1 = cos(pi/(N-1)*(0:(N-1))'*(0:(N-1)));
DCT1(:,[1,end]) = .5*DCT1(:,[1,end]);
y = DCT1 * x;
pass(5) = norm( y - chebfun.dct(x, 1) ) < tol; 

% DCT-II
DCT2 = cos(pi/N*(0:N-1)'*((0:N-1)+1/2));
y = DCT2 * x; 
pass(6) = norm( y - chebfun.dct(x, 2) ) < tol; 

% DCT-III
DCT3 = cos(pi/N*((0:N-1)+1/2)'*(0:N-1));
DCT3(:,1) = .5*DCT3(:,1);
y = DCT3 * x; 
pass(7) = norm( y - chebfun.dct(x, 3) ) < tol; 

% DCT-IV
DCT4 = cos(pi/N*((0:N-1)+1/2)'*((0:N-1)+1/2));
y = DCT4 * x; 
pass(8) = norm( y - chebfun.dct(x, 4) ) < 10*tol; 

% IDCT-I 
y = DCT1 \ x;
pass(9) = norm( y - chebfun.idct(x, 1) ) < tol; 

% IDCT-II
y = DCT2 \ x; 
pass(10) = norm( y - chebfun.idct(x, 2) ) < tol; 

% IDCT-III
y = DCT3 \ x; 
pass(11) = norm( y - chebfun.idct(x, 3) ) < tol; 

% IDCT-IV
y = DCT4 \ x; 
pass(12) = norm( y - chebfun.idct(x, 4) ) < tol; 

end
