function d = diffbarytrig(zz, zj, fj, wj, varargin)
%DIFFBARYTRIG  Computes derivative of trigonometric rational function in
%barycentric form.
%   D = DIFFBARYTRIG(ZZ, ZJ, FJ, WJ) returns D (vector of floats), the values of the
%   derivative of a odd trigonometric barycentric rational function with support points ZJ,
%   function values FJ, and barycentric weights WJ evaluated at the points ZZ.
%
%   D = DIFFBARYTRIG(ZZ, ZJ, FJ, WJ, N) computes the N-th derivative.
%  
%   D = DIFFBARYTRIG(ZZ, ZJ, FJ, WJ, FORM) computes derivative of a
%   barycentric rational function of type FORM, where FORM can be 'odd' or
%   'even'. The default value of FORM is 'odd'.
%
% See also AAATRIG, REVALTRIG.

% Copyright 2018 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

form = 'odd'; % Default value for form
N = 1; % Default value for derivative

while (~isempty(varargin))
    isinteger(varargin{1})
    if  strncmpi(varargin{1},'even',4)
          form = 'even';
          varargin(1) = [];
    elseif strncmpi(varargin{1},'odd',3)
          form = 'odd';
          varargin(1) = []; 
    elseif isfloat(varargin{1})
          N = varargin{1};
          varargin(1) = [];
          if N==0; error('Use revaltrig.m for function evaluation'); end
    end
end

m = numel(zj);
npts = numel(zz);

zv = zz(:);                            % Vectorize zz if necessary
zvp= zv - 2*pi*floor(real(zv/(2*pi))); % Project evaluation points onto first period window

% Define basis functions
if strcmp(form,'even')
    cst = @(z) cot(z);
elseif strcmp(form,'odd')
    cst = @(z) csc(z);
end    

CC = cst(bsxfun(@minus, zvp, zj.')/2); % Cauchy matrix
rpDen = CC*wj; % Denominator of derivative
rn = (CC*(wj.*fj))./(CC*wj);             % vector of values

rp = zeros(npts,N+1);
D = zeros(m,m,N+1);
D(:,:,1) = eye(m);
rp(:,1) = rn;

% Compute the N derivatives
for p = 1:N
    % Derivative away from support points
    DR = zeros(npts,m,p);
    q = permute(0:(p-1),[1,3,2]);
    DR(:,:,q+1) = factorial(p)./(factorial(q).*factorial(p-q))...
                .*2.^(q-p).*diffCst((zvp - zj.')/2,p-q,form).*...
                (heaviside(.5-q).*fj.' - permute(rp(:,q+1),[1,3,2])).*wj.';
    rp(:,p+1) = sum(sum(DR,2),3)./rpDen;    

    q = permute(1:p,[1,3,2]);
    
    % Now compute the differentiation matrix
    % Compute the terms in the first sum of the formula
    firstSum =  wj.'./wj.*sum(bsxfun(@times,eye(m),D(:,:,flip(1:p))),2)...
                .*2.^(-q).*diffCstInv((zj - zj)/2,q,form);
    % Compute the terms in the second sum in the formula        
    secondSum =  D(:,:,flip(1:p)).*2.^(-q).*diffCstInv((zj - zj.')/2,q,form);
    % Combine the sums with the coefficients
    D1 = cst((zj-zj.')/2).*sum(factorial(p)./(factorial(q).*factorial(p-q))...
               .*(firstSum - secondSum),3);
    % Replace NaN terms with correct terms from formula       
    D1(1:m+1:m^2) = 0;
    D1 = D1 - diag(sum(D1(:,:),2));
    D(:,:,p+1) = D1;
  
end

ii = find(isnan(rp(:,N+1)));
DZJ = D(:,:,end)*fj; % derivative at zj

% Correct NaN derivatives, including derivative at support points
for jj = 1:length(ii)
    if ( isnan(zvp(ii(jj))) || ~any(zvp(ii(jj)) == zj) )
        % r(NaN) = NaN is fine.
        % The second case may happen if r(zv(ii)) = 0/0 at some point.
    else
        % Clean up values NaN = inf/inf at support points.
        % Find the corresponding node and set entry to correct value:
        rp(ii(jj),N+1) = DZJ(zvp(ii(jj)) == zj);
    end
end

% Deal with input inf
rp(isinf(real(zvp/1i)),N+1) = 0;

% Reshape to input format:
d = reshape(rp(:,N+1), size(zz));

end %

function d = diffCot(t,n) % Returns the n-th derivative of Cot evaluated at t
% using derivative polynomials. The t can take arbitrary dimension.

x = tan(t + pi/2);
Psz = [size(x), max(n)+1];
nDim = numel(Psz); % Define the final dimension for summation
P = zeros(Psz); % Initialize derivativ

% Define indices for summations
inds1 = repmat({':'},nDim,1);
inds0 = repmat({1},nDim,1);
inds2 = inds1;
inds3 = inds1;

inds1{nDim}=1;
P(inds1{:}) = -x;

for k = 0:(max(n)-1)
    inds0{numel(Psz)}=k+1;
    l = reshape((0:k)',inds0{:});
    inds1{nDim}=k+2;
    inds2{nDim}=l+1;
    inds3{nDim}=1+k-l;

    P(inds1{:}) = -sum(factorial(k)./(factorial(l).*factorial(k-l)).*P(inds2{:}).*P(inds3{:}),numel(Psz))...
              - heaviside(eps-k);
end

inds1{nDim}=n+1;
d = P(inds1{:});

end

function d= diffCsc(t,n)

d = (0.5).^n.*diffCot(t/2,n) - diffCot(t,n);

end

function d= diffSin(t,n)

d = sin(t+n*pi/2);

end

function d= diffTan(t,n)

d =  -diffCot(t-pi/2,n);

end

function d = diffCstInv(t,n,form)

if strcmp(form,'even')
    d = diffTan(t,n);
elseif strcmp(form,'odd')
    d = diffSin(t,n);
end

end

function d = diffCst(t,n,form)

if strcmp(form,'even')
    d = diffCot(t,n);
elseif strcmp(form,'odd')
    d = diffCsc(t,n);
end

end