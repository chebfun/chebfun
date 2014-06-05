function C = mtimes(A, B)
%MTIMES   Composition of LINOPs.

C = mtimes@chebmatrix(A,B);

if ( isa(A,'linop') && isa(B,'linop') )
    C = linop(C);
    C.hasGivenJumpsAt = union(A.hasGivenJumpsAt,B.hasGivenJumpsAt);
end

end
