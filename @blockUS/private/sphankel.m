function H = sphankel(r)
% SPHANKEL(R) this forms a sparse hankel matrix by forming it as an upside-
% down toeplitz matrix. 

%Hankel is an upside-down toeplitz matrix. 
r = flipud(r(:));   %ensure column vector. 
H = fliplr(triu(sptoeplitz(r,r)));
end