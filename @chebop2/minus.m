function N = minus(N1, N2)
%MINUS minus for chebop2 objects

N = plus( N1, uminus(N2) );

end