function [pol, res, zer] = prztrig(zj, fj, wj, form)
%PRZTRIG   Computes poles, residues, and zeros of a periodic rational function
%          in trigonometric barycentric form.
%   [POL, RES, ZER] = PRZTRIG(ZJ, FJ, WJ) returns vectors of poles POL,
%   residues RES, and zeros ZER of the periodic rational function defined by
%   support points ZJ, function values FJ, and barycentric weights WJ.
%
% See also AAATRIG, PRZ, REVALTRIG.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

m = length(wj);

if strcmp(form,'even') % Even (cot) approximant
d = (zj ~= pi); % Special case with support point at half of period
zjp = tan(zj(d)/2); % Coordinate transformation.
wjp = wj(d).*(1+zjp.^2);
if abs(zjp)>1e5; warning('The transformation for finding the poles and zeros may be ill-conditioned.'); end
    
% Compute poles via generalized eigenvalue problem:
B = eye(m+1);
B(1,1) = 0;
cd = sum(zjp.*wj(d)); % Constant term due to transformation

if all(d) % When none of the support points are pi  
Ep = [cd wjp.'; ones(m, 1) diag(zjp)];    
else % When one of the support points is pi
Ep = [-wj(~d) cd wjp.'; 
      1 zeros(1,m);
      zeros(m-1,1) ones(m-1, 1) diag(zjp)];
Ep(2,2) = 0;
end
    
% Find poles as eigenvalues
[vpolp,polp] = eig(Ep, B);
polp = diag(polp);
polp = polp(~isinf(polp)); % Remove spurious pole(s)
pol = 2*atan(polp); % Invert transformation

% Compute zeros via generalized eigenvalue problem:
cn = sum(fj(d).*zjp.*wj(d)); % Constant term due to transformation

if all(d) % When none of the support points are pi
Ez = [cn fj(d).'.*wjp.'; ones(m, 1) diag(zjp)];    
else % When one of the support points is pi
Ez = [-fj(~d)*wj(~d)  cn  fj(d).'.*wjp.'; 
    1 zeros(1,m);
    zeros(m-1,1) ones(m-1, 1) diag(zjp)];
Ez(2,2) = 0;
end

% Find zeros as eigenvalues
[vzerp,zerp] = eig(Ez, B);
zerp = diag(zerp);
zerp = zerp(~isinf(zerp)); % Remove spurious zero(s)
zer = 2*atan(zerp); % Reverse transformation

% Account for poles and zeros at infinity
zer(abs(zerp-1i)<1e-10) =  1i*Inf;
zer(abs(zerp+1i)<1e-10) = -1i*Inf;
pol(abs(polp-1i)<1e-10) =  1i*Inf;
pol(abs(polp+1i)<1e-10) = -1i*Inf;

else % Odd (csc) approximant

zjp = exp(1i*zj); % Coordinate transformation.
wjp = wj.*exp(1i*zj/2);

% Compute poles via generalized eigenvalue problem:
B = eye(m+1);
B(1,1) = 0;

Ep = [0 wjp.'; ones(m, 1) diag(zjp)];    
% Find poles as eigenvalues
[vpolp,polp] = eig(Ep, B);
polp = diag(polp);
polp = polp(~isinf(polp)); % Remove spurious pole(s)
pol = -1i*log(polp); % Invert transformation

Ez = [0 fj.'.*wjp.'; ones(m, 1) diag(zjp)];    
% Find zeros as eigenvalues
[vzerp,zerp] = eig(Ez, B);
zerp = diag(zerp);
zerp = zerp(~isinf(zerp)); % Remove spurious zero(s)
zer = -1i*log(zerp); % Reverse transformation

% Account for poles and zeros at infinity
zer(abs(zerp)<1e-10) =  1i*Inf;
zer(abs(zerp)>1e10) = -1i*Inf;
pol(abs(polp)<1e-10) =  1i*Inf;
pol(abs(polp)>1e10) = -1i*Inf;

end

% Remove pairs poles and zeros that coincide at +1i*Infinity
zerPInf = find(zer == 1i*Inf);
polPInf = find(pol == 1i*Inf);
numPInf = min(length(zerPInf),length(polPInf)); % Number of coincident poles and zeros
zer(zerPInf(1:numPInf)) = [];
pol(polPInf(1:numPInf)) = [];

% Now do the same this for -1i*Infinity
zerMInf = find(zer == -1i*Inf);
polMInf = find(pol == -1i*Inf);
numMInf = min(length(zerMInf),length(polMInf)); % Number of coincident poles and zeros
zer(zerMInf(1:numMInf)) = [];
pol(polMInf(1:numMInf)) = [];

% Project poles and zeros into first period window
pol = pol - 2*pi*floor(real(pol/(2*pi)));
zer = zer - 2*pi*floor(real(zer/(2*pi)));    

% Now calculate residues
if ~isempty(pol)
    % Compute residues via formula for res of quotient of analytic functions:
    if strcmp(form,'even') % Even (cot) approximant
        N = @(t) (cot(bsxfun(@minus,t,zj.')/2)) * (fj.*wj);
        D = @(t) (cot(pi/P*bsxfun(@minus,t,zj.')/2)) * wj;
        Ddiff = @(t) -1/2*csc(bsxfun(@minus,t,zj.')/2).^2 * wj;
    elseif strcmp(form,'odd') % Odd (csc) approximant
        % Compute residues via formula for res of quotient of analytic functions:
        N = @(t) (csc(bsxfun(@minus,t,zj.')/2)) * (fj.*wj);
        D = @(t) (csc(pi/P*bsxfun(@minus,t,zj.')/2)) * wj;
        Ddiff = @(t) -1/2*csc(bsxfun(@minus,t,zj.')/2).*cot(bsxfun(@minus,t,zj.')/2) * wj;
    end
    res = N(pol)./Ddiff(pol);
else
    res = [];
end

end % End of PRZTRIG().
