function F = uminus( F )
% + Unary minus of a chebfun2v.

components = F.components; 
for j = 1:F.nComponents
    components(j) = -components(j);
end
F.components = components; 
end