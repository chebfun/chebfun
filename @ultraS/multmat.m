function M = multmat(n, f, lambda)
% MULTMAT  multiplication matrices for ultraS
%
%  M = MULTMAT(N, F, LAMBDA) forms the nxn multiplication matrix
%  representing the multiplication of F in the C^{(LAMBDA)} basis.
% 
%  M = MULTMAT(N, F, LAMBDA) also works when F is a vector of Chebyshev
%  coefficients.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( isa(f, 'chebfun') ) 
    % Get Chebyshev T coefficients
    a = get(f, 'coeffs');
elseif ( isa(f, 'bndfun') ) 
    % Get Chebyshev T coefficients
    a = get(f, 'coeffs');
elseif ( isa( f, 'double') )
    a = f; 
else
    error('ULTRAS:ARGIN:TYPE', 'Unrecognised 2nd argument.');
end

% Multiplying by a scalar is easy.
if ( numel(a) == 1 )
    M = a*speye(n);
    return
end

% Prolong or truncate coefficients
if ( numel(a) < n )
    a = [a ; zeros(n - numel(a), 1)];   % Prolong
else
    a = a(1:n);                         % Truncate.
end

if ( lambda == 0 )          % Multiplication in Chebyshev T coefficients.
    a = a/2;  % just to make formula easier.
    M = ultraS.sptoeplitz([2*a(1);a(2:end)], [2*a(1);a(2:end)]);
    H = ultraS.sphankel(a(2:end));
    % TODO: What do the variables sub1 and sub2 represent?
    sub1 = 2:length(a); 
    sub2 = 1:length(a)-1;
    M(sub1, sub2) = M(sub1, sub2) + H;
    
elseif ( lambda == 1 )      % Multiplication in Chebyshev U coefficients.
    M = ultraS.sptoeplitz([2*a(1);a(2:end)], [2*a(1);a(2:end)])/2;
    sub = 1:length(a) - 2;
    M(sub, sub) = M(sub, sub) - ultraS.sphankel(a(3:end)/2);
    
else
    % Want the C^{lam}C^{lam} Cheb Multiplication matrix.
   
    % Convert ChebT of a to ChebC^{lam}
    a = ultraS.convertmat(n, 0, lambda - 1) * a;
    
    m = 2*n; 
    M0 = speye(m);
    
    d1 = [1 (2*lambda : 2*lambda + m - 2)]./ ...
        [1 (2*((lambda+1) : lambda + m - 1))];
    d2 = (1:m)./(2*(lambda:lambda + m - 1));
    B = [d2' zeros(m, 1) d1'];
    Mx = spdiags(B,[-1 0 1], m, m);
    M1 = 2*lambda*Mx;
    
    % Construct the multiplication operator by a three-term recurrence: 
    M = a(1)*M0;
    M = M + a(2)*M1;
    for nn = 1:length(a) - 2
        M2 = 2*(nn + lambda)/(nn + 1)*Mx*M1 - (nn + 2*lambda - 1)/(nn + 1)*M0;
        M = M + a(nn + 2)*M2;
        M0 = M1;
        M1 = M2;
        if ( abs(a(nn + 3:end)) < eps ), break, end
    end
    M = M(1:n, 1:n); 
end

end
