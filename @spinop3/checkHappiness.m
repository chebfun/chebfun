function ishappy = checkHappiness(S, c, pref)
%CHECKHAPPINESS   Check if the solution is resolved in space when using SPIN3.
%   ISHAPPY = CHECKHAPPINESS(S, C, PREF) checks if the solution C of the PDE 
%   defined by the SPINOP3 S is resolved in space relative to PREF.ERRTOL.
%   ISHAPPY is 1 if it is resolved, 0 otherwise.
%
% See also SPIN3.

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
    cCheb = ifftshift(c{1}(idx:idx+N-1,setToCheck,setToCheck));
    if ( mod(N,2) == 0 )
        cCheb = [cCheb(N); cCheb(N-1:-1:N/2+1) + cCheb(1:N/2-1); cCheb(N/2)];
    else
        cCheb = [cCheb(N:-1:(N+1)/2+1) + cCheb(1:(N+1)/2-1); cCheb((N+1)/2)];
    end
    
    % Trick used in TRIGTECH/STANDARCHECK to get decaying coefficients:
    cCheb = flipud(cCheb);
    cCheb = [cCheb(1,1); kron(cCheb(2:end,1), [1; 1])];
    
    % Use CHEBFUN STANDARDCHOP command:
    cutoff = standardChop(cCheb, errTol);
    ishappy(k) = ( cutoff < N );
    
    % Second set of coefficients.
    % IFFTSHIFT to use CHEBFUN indexing of wavenumbers:
    setToCheck = N;
    cCheb = ifftshift(c{1}(idx:idx+N-1,setToCheck,setToCheck));
    if ( mod(N,2) == 0 )
        cCheb = [cCheb(N); cCheb(N-1:-1:N/2+1) + cCheb(1:N/2-1); cCheb(N/2)];
    else
        cCheb = [cCheb(N:-1:(N+1)/2+1) + cCheb(1:(N+1)/2-1); cCheb((N+1)/2)];
    end
    
    % Trick used in TRIGTECH/STANDARCHECK to get decaying coefficients:
    cCheb = flipud(cCheb);
    cCheb = [cCheb(1,1); kron(cCheb(2:end,1), [1; 1])];
    
    % Use CHEBFUN STANDARDCHOP command:
    cutoff = standardChop(cCheb, errTol);
    ishappy(k) = ishappy(k) && ( cutoff < N );
    
end
ishappy = all(ishappy == 1);

end