function F = ctranspose( F )
% '   Conjugate transpose of a CHEBFUN3V

% Transpose and then conjugate: 
F = transpose( F ); 
F = conj( F ); 

end
