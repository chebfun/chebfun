function pass = test_dct(pref)

if ( nargin == 0 )
    pref = chebfunpref();
end

% TODO: include more extensive tests.

n = 5;
seedRNG(42)
R = rand(n);

tol = n*eps;

for k = 1:8
    pass(k) = norm(chebfun.idct(chebfun.dct(R,1),1) - R) < tol;
end

N = 13; 
x = rand(N, 2); 
tol = 50*eps;

% DCT-I
DCT1 = cos(pi/(N-1)*(0:(N-1))'*(0:(N-1)));
DCT1(:,[1,end]) = .5*DCT1(:,[1,end]);
y = DCT1 * x;
pass(9) = norm( y - chebfun.dct(x, 1) ) < tol; 

% DCT-II
DCT2 = cos(pi/N*(0:N-1)'*((0:N-1)+1/2));
y = DCT2 * x; 
pass(10) = norm( y - chebfun.dct(x, 2) ) < tol; 

% DCT-III
DCT3 = cos(pi/N*((0:N-1)+1/2)'*(0:N-1));
DCT3(:,1) = .5*DCT3(:,1);
y = DCT3 * x; 
pass(11) = norm( y - chebfun.dct(x, 3) ) < tol; 

% DCT-IV
DCT4 = cos(pi/N*((0:N-1)+1/2)'*((0:N-1)+1/2));
y = DCT4 * x; 
pass(12) = norm( y - chebfun.dct(x, 4) ) < 10*tol; 

% DCT-V 
DCT5 = cos(pi*(0:N-1)'*(0:N-1)/(N-1/2));
y = DCT5 * x; 
pass(13) = norm( y - chebfun.dct(x, 5) ) < 10*tol; 

% DCT-VI 
DCT6 = cos(pi*(0:N-1)'*((0:N-1)+1/2)/(N-1/2));
y = DCT6 * x; 
pass(14) = norm( y - chebfun.dct(x, 6) ) < 10*tol;

% DCT-VII 
DCT7 = cos(pi*((0:N-1)+1/2)'*(0:N-1)/(N-1/2));
y = DCT7 * x; 
pass(15) = norm( y - chebfun.dct(x, 7) ) < 10*tol; 

% DCT-VIII
DCT8 = cos(pi*((0:N-1)+1/2)'*((0:N-1)+1/2)/(N+1/2));
y = DCT8 * x; 
pass(16) = norm( y - chebfun.dct(x, 8) ) < 10*tol; 

% IDCT-I 
y = DCT1 \ x;
pass(17) = norm( y - chebfun.idct(x, 1) ) < tol; 

% IDCT-II
y = DCT2 \ x; 
pass(18) = norm( y - chebfun.idct(x, 2) ) < tol; 

% IDCT-III
y = DCT3 \ x; 
pass(19) = norm( y - chebfun.idct(x, 3) ) < tol; 

% IDCT-IV
y = DCT4 \ x; 
pass(20) = norm( y - chebfun.idct(x, 4) ) < tol; 

% IDCT-V
y = DCT5 \ x;
pass(21) = norm( y - chebfun.idct(x, 5) ) < tol; 

% IDCT-VI
y = DCT6 \ x; 
pass(22) = norm( y - chebfun.idct(x, 6) ) < tol; 

% IDCT-VII
y = DCT7 \ x; 
pass(23) = norm( y - chebfun.idct(x, 7) ) < tol; 

% IDCT-VIII
y = DCT8 \ x; 
pass(24) = norm( y - chebfun.idct(x, 8) ) < tol; 


% Check complex inputs 
c = rand(10,1) + 1i*rand(10,1); 
for j = 1:8
    v1 = chebfun.dct(c, j);
    v2 = chebfun.dct(real(c),j) + 1i*chebfun.dct(imag(c),j);
    pass(24+j) = norm(v1 - v2) < tol; 
end

end
