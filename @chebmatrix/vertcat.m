function C = vertcat(A,B)
C = chebmatrix( vertcat(A.blocks,B.blocks) );
end
