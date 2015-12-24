function v = norm( F )
%NORM   Frobenius norm of a SPHEREFUNV.
%       V = sqrt(norm(F1).^2 + norm(F2).^2 + norm(F3).^2) .

% Empty check: 
if ( isempty( F ) ) 
    v = []; 
    return
end

v = 0; 
for jj = 1:3 
    v = v + sum( norm( F.components{jj}, 2 ).^2 );
end
v = sqrt(v); 

end
