function C = minus(A, B)
%MINUS   Difference of LINOPs.

C = linop( minus@chebmatrix(A,B) );
C.hasGivenJumpsAt = union(A.hasGivenJumpsAt,B.hasGivenJumpsAt);

end
