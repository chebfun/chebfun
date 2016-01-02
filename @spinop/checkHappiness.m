function ishappy = checkHappiness(S, c, pref)
%CHECKHAPPINESS   Check if the solution is resolved in space for SPINOP.
%   CHECKHAPPINESS 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
nVars = S.numVars;
N = size(c{1}, 1)/nVars;
errTol = pref.errTol;

% Loop over the variables:
ishappy = zeros(nVars, 1);
for k = 1:nVars
    
    idx = (k-1)*N + 1;
    
    % IFFTSHIFT to use CHEBFUN indexing of wavenumbers:
    cShift = ifftshift(c{1}(idx:idx+N-1));
    if ( mod(N,2) == 0 )
        cShift = [cShift(N); cShift(N-1:-1:N/2+1) + cShift(1:N/2-1); ...
            cShift(N/2)];
    else
        cShift = [cShift(N:-1:(N+1)/2+1,:) + cShift(1:(N+1)/2-1,:); ...
            cShift((N+1)/2,:)];
    end
    
    % Trick used in TRIGTECH/STANDARCHECK to get decaying coefficients:
    cShift = flipud(cShift);
    cShift = [cShift(1,1); kron(cShift(2:end,1), [1; 1])];
    
    % Use CHEBFUN STANDARDCHOP command:
    cutoff = standardChop(cShift, errTol);
    ishappy(k) = ( cutoff < N );
    
end
ishappy = all(ishappy == 1);

end