function ishappy = checkHappiness(S, c, pref)
%CHECKHAPPINESS   Check if the solution is resolved in space for SPINIOP3.
%   CHECKHAPPINESS 

% Copyright 2016 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

% Set-up:
nVars = S.numVars;
N = size(c{1}, 1)/nVars;
errTol = pref.errTol;

% Loop over the variables.
% Remark: We only check two sets of coefficients per variable.
ishappy = zeros(nVars, 1);
for k = 1:nVars
    
    idx = (k-1)*N + 1;
    
    % First set of coefficients.
    % IFFTSHIFT to use CHEBFUN indexing of wavenumbers:
    setToCheck = 1;
    cShift = ifftshift(c{1}(idx:idx+N-1,setToCheck,setToCheck));
    if ( mod(N,2) == 0 )
        cShift = [cShift(N); cShift(N-1:-1:N/2+1) + ...
            cShift(1:N/2-1); cShift(N/2)];
    else
        cShift = [cShift(N:-1:(N+1)/2+1,:) + ...
            cShift(1:(N+1)/2-1,:); cShift((N+1)/2,:)];
    end
    
    % Trick used in TRIGTECH/STANDARCHECK to get decaying coefficients:
    cShift = flipud(cShift);
    cShift = [cShift(1,1); kron(cShift(2:end,1), [1; 1])];
    
    % Use CHEBFUN STANDARDCHOP command:
    cutoff = standardChop(cShift, errTol);
    ishappy(k) = ( cutoff < N );
    
    % Second set of coefficients.
    % IFFTSHIFT to use CHEBFUN indexing of wavenumbers:
    setToCheck = N;
    cShift = ifftshift(c{1}(idx:idx+N-1,setToCheck,setToCheck));
    if ( mod(N,2) == 0 )
        cShift = [cShift(N); cShift(N-1:-1:N/2+1) + ...
            cShift(1:N/2-1); cShift(N/2)];
    else
        cShift = [cShift(N:-1:(N+1)/2+1,:) + ...
            cShift(1:(N+1)/2-1,:); cShift((N+1)/2,:)];
    end
    
    % Trick used in TRIGTECH/STANDARCHECK to get decaying coefficients:
    cShift = flipud(cShift);
    cShift = [cShift(1,1); kron(cShift(2:end,1), [1; 1])];
    
    % Use CHEBFUN STANDARDCHOP command:
    cutoff = standardChop(cShift, errTol);
    ishappy(k) = ishappy(k) && ( cutoff < N );
    
end
ishappy = all(ishappy == 1);

end