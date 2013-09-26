function C = horzcat(A,B)
if isa(A,'chebmatrix')
    b1 = A.blocks;
else
    b1 = {A};
end
if isa(B,'chebmatrix')
    b2 = B.blocks;
else
    b2 = {B};
end
C = chebmatrix( horzcat(b1,b2) );
end
