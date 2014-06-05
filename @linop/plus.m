function C = plus(A, B)
%PLUS   Addition of LINOPs.

C = linop( plus@chebmatrix(A,B) );
C.hasGivenJumpsAt = union(A.hasGivenJumpsAt,B.hasGivenJumpsAt);

end
