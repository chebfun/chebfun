function cf3f = constructor(cf3f, f, varargin)

% Parse the input
pref             = chebfunpref();
dom = [-1, 1, -1, 1, -1, 1];
for k = 1:length(varargin)
    if strcmpi(varargin{k}, 'eps')
        pref.cheb3Prefs.chebfun3eps = varargin{k+1};
    end
    if ( numel(varargin{k}) == 6 ) % domain is specified.
        dom = varargin{k};
    end
end

% Set Chebfun parameters:
tech             = pref.tech();
prefStruct       = pref.cheb3Prefs;
tpref            = tech.techPref;
grid             = tpref.minSamples;
maxSample        = tpref.maxLength;     % max polynomialDeg
maxSamplePhase1  = 363;                 % max coarseResolution (not implemented)
maxRank          = prefStruct.maxRank;  % max rank (not implemented)
pseudoLevel      = prefStruct.chebfun3eps;
passSampleTest   = prefStruct.sampleTest;
maxRestarts      = 10;


% Check if fast vectorized evaluations are possible
try
    A = f(1:2,1:2,1:2);
    vectorize = 1;
    if isscalar(A)
        f = @(x,y,z) f(x,y,z) + 0*x + 0*y + 0*z; %ensure vectorization
    end
catch
    vectorize = 0;
end

% Initialize
n                = [grid, grid, grid];  %coarseResolution
m                = n;                   %fineResolution
r                = [6,6,6];             %rank
tol              = pseudoLevel;
chebX            = @(i,n) dom(1) + ((-cos((i-1).*pi/(n-1))) + 1)*(dom(2)-dom(1))/2 ;
chebY            = @(i,n) dom(3) + ((-cos((i-1).*pi/(n-1))) + 1)*(dom(4)-dom(3))/2 ;
chebZ            = @(i,n) dom(5) + ((-cos((i-1).*pi/(n-1))) + 1)*(dom(6)-dom(5))/2 ;
reffun           = @(n) floor(sqrt(2)^(floor(2*log2(n)) + 1)) + 1;
restarts         = 0;
cf3f.domain      = dom;

%% Main Loop
happy = 0;
while ~happy
    
    %% Phase 1
    happyPhase1 = 0;
    while ~happyPhase1
        J = initializeIndexRandomly(r(2), n(2));
        K = initializeIndexRandomly(r(3), n(3));
        
        % Handle to evaluate tensor entries of T_c
        ff = @(i,j,k) f(chebX(i,n(1)),chebY(j,n(2)),chebZ(k,n(3)));
        
        for iterations = 1:2
            happyPhase1 = 1;
            
            % ACA 1
            T1 = evalTensor(1:n(1),J,K,ff,vectorize);
            T1 = reshape(T1,n(1),r(2)*r(3));
            [~, tol] = getTol(T1, pseudoLevel, tol,dom(2)-dom(1));
            [Uc, ~, ~, I,I2] = ACA(T1, tol, n(1));
            r(1) = size(I,2);
            JT1 = J;
            KT1 = K;
            
            % ACA 2
            T2 = evalTensor(I,1:n(2),K,ff,vectorize);
            T2 = reshape(permute(T2,[2,1,3]),n(2),r(1)*r(3));
            [~, tol] = getTol(T2, pseudoLevel, tol,dom(4)-dom(3));
            [Vc, ~, ~, J, J2] = ACA(T2, tol, n(2));
            r(2) = size(J,2);
            KT2 = K;
            IT2 = I;
            
            % ACA 3
            T3 = evalTensor(I,J,1:n(3), ff,vectorize);
            T3 = reshape(permute(T3,[3,1,2]),n(3),r(1)*r(2));
            [reltol, tol] = getTol(T3, pseudoLevel, tol,dom(6)-dom(5));
            [Wc, ~, ~, K, K2] = ACA(T3, tol, n(3));
            r(3) = size(K,2);
            IT3 = I;
            JT3 = J;
            
            % Refine n
            breakFlag = 0;
            while r(1)*2*sqrt(2) > n(1)
                n(1) = reffun(n(1));
                breakFlag = 1;
            end
            while r(2)*2*sqrt(2) > n(2)
                n(2) = reffun(n(2));
                for i = 1:3
                    while r(i)*2*sqrt(2) > n(i)
                        n(i) = reffun(n(i));
                    end
                end
                breakFlag = 1;
            end
            while r(3)*2*sqrt(2) > n(3)
                n(3) = reffun(n(3));
                breakFlag = 1;
            end
            
            % Reinitialize after refinement
            if breakFlag == 1
                happyPhase1 = 0;
                break
                % Proceed to refinement if r < 2
            elseif min(r) < 2
                break
            end
        end
    end
    
    % Catch the rank zero function:
    if size(I,1) == 0 || size(J,1) == 0 || size(K,1) == 0
        cf3f.cols = chebfun(zeros([n(1),1]), [dom(1),dom(2)], pref);
        cf3f.rows = chebfun(zeros([n(2),1]), [dom(3),dom(4)], pref);
        cf3f.tubes = chebfun(zeros([n(3),1]), [dom(5),dom(6)], pref);
        cf3f.core = 0;
    else
        %% Refine
        pref.chebfuneps = reltol;
        m = n;
        
        % Check if further refinement is necessary
        %U
        Uf = Uc;
        fiberData.hscale = norm(dom(1:2), inf);
        ct2 = createCT2(Uf,fiberData);
        resolvedU = happinessCheck(ct2, [], ct2.coeffs, [], pref);
        if ~resolvedU
            m(1) = 2*m(1)-1;
        end
        
        %V
        Vf = Vc;
        fiberData.hscale = norm(dom(3:4), inf);
        ct2 = createCT2(Vf, fiberData);
        resolvedV = happinessCheck(ct2, [], ct2.coeffs, [], pref);
        if ~resolvedV
            m(2) = 2*m(2)-1;
        end
        
        % W
        Wf = Wc;
        fiberData.hscale = norm(dom(5:6), inf);
        ct2 = createCT2(Wf,fiberData);
        resolvedW = happinessCheck(ct2, [], ct2.coeffs, [], pref);
        if ~resolvedW
            m(3) = 2*m(3)-1;
        end
        
        % Add function evaluations and check again
        while ~resolvedU || ~resolvedV || ~resolvedW
            % function handle to evaluate T_f
            ff = @(i,j,k) f(chebX(i,m(1)),chebY(j,m(2)),chebZ(k,m(3)));
            
            % Map the indices from T_c to T_f
            refFactor = [0, 0, 0];
            iter = 0;
            while min(refFactor) == 0
                iter = iter +1;
                if n(1)*iter-(iter-1) == m(1)
                    refFactor(1) = iter;
                end
                if n(2)*iter-(iter-1) == m(2)
                    refFactor(2) = iter;
                end
                if n(3)*iter-(iter-1) == m(3)
                    refFactor(3) = iter;
                end
            end
            ref = @(i, r) r*i-(r-1);
            
            % U
            Jr = ref(JT1, refFactor(2));
            Kr = ref(KT1, refFactor(3));
            if ~resolvedU
                Uold = Uf;
                [Ij,Ik] = ind2sub([size(Jr,2),size(Kr,2)],I2);
                if vectorize == 1
                    nn = [floor(m(1)/2),size(I2,2)];
                    x = zeros([nn(1),1]);
                    x(:,1) = 2:2:m(1);
                    X = repmat(x,1,nn(2));
                    y = zeros([1,nn(2)]);
                    y(1,:) = Jr(Ij);
                    Y = repmat(y,nn(1),1);
                    z = zeros([1,nn(2)]);
                    z(1,:) = Kr(Ik);
                    Z = repmat(z,nn(1),1);
                    Uf(2:2:m(1),1:nn(2)) = ff(X,Y,Z);
                    Uf(1:2:m(1),1:nn(2)) = Uold;
                else
                    Uf = zeros([m(1),size(I2,2)]);
                    for i = 1:m(1)
                        for j = 1:size(I2,2)
                            
                            if mod(i,2) == 0
                                Uf(i,j) = ff(i, Jr(Ij(j)), Kr(Ik(j)));
                            else
                                Uf(i,j) = Uold((i+1)/2,j);
                            end
                        end
                    end
                end
                fiberData.hscale = norm(dom(1:2), inf);
                ct2 = createCT2(Uf,fiberData);
                resolvedU = happinessCheck(ct2, [], ct2.coeffs, [], pref);
                if ~resolvedU
                    m(1) = 2*m(1)-1;
                end
            end
            
            % V
            Ir = ref(IT2, refFactor(1));
            Kr = ref(KT2, refFactor(3));
            if ~resolvedV
                Vold = Vf;
                [Ji,Jk] = ind2sub([size(Ir,2),size(Kr,2)],J2);
                if vectorize == 1
                    nn = [floor(m(2)/2),size(J2,2)];
                    y = zeros([nn(1),1]);
                    y(:,1) = 2:2:m(2);
                    Y = repmat(y,1,nn(2));
                    x = zeros([1,nn(2)]);
                    x(1,:) = Ir(Ji);
                    X = repmat(x,nn(1),1);
                    z = zeros([1,nn(2)]);
                    z(1,:) = Kr(Jk);
                    Z = repmat(z,nn(1),1);
                    Vf(2:2:m(2),1:nn(2)) = ff(X,Y,Z);
                    Vf(1:2:m(2),1:nn(2)) = Vold;
                else
                    Vf = zeros([m(2),size(J2,2)]);
                    for i = 1:m(2)
                        for j = 1:size(J2,2)
                            if mod(i,2) == 0
                                Vf(i,j) = ff(Ir(Ji(j)), i, Kr(Jk(j)));
                            else
                                Vf(i,j) = Vold((i+1)/2,j);
                            end
                        end
                    end
                end
                fiberData.hscale = norm(dom(3:4), inf);
                ct2 = createCT2(Vf,fiberData);
                resolvedV = happinessCheck(ct2, [], ct2.coeffs, [], pref);
                if ~resolvedV
                    m(2) = 2*m(2)-1;
                end
            end
            
            % W
            Ir = ref(IT3, refFactor(1));
            Jr = ref(JT3, refFactor(2));
            if ~resolvedW
                Wold = Wf;
                [Ki,Kj] = ind2sub([size(Ir,2),size(Jr,2)],K2);
                if vectorize == 1
                    nn = [floor(m(3)/2),size(K2,2)];
                    z = zeros([nn(1),1]);
                    z(:,1) = 2:2:m(3);
                    Z = repmat(z,1,nn(2));
                    x = zeros([1,nn(2)]);
                    x(1,:) = Ir(Ki);
                    X = repmat(x,nn(1),1);
                    y = zeros([1,nn(2)]);
                    y(1,:) = Jr(Kj);
                    Y = repmat(y,nn(1),1);
                    Wf(2:2:m(3),1:nn(2)) = ff(X,Y,Z);
                    Wf(1:2:m(3),1:nn(2)) = Wold;
                else
                    Wf = zeros([m(3),size(K2,2)]);
                    for i = 1:m(3)
                        for j = 1:size(K2,2)
                            if mod(i,2) == 0
                                Wf(i,j) = ff(Ir(Ki(j)), Jr(Kj(j)), i);
                            else
                                Wf(i,j) = Wold((i+1)/2,j);
                            end
                        end
                    end
                end
                fiberData.hscale = norm(dom(5:6), inf);
                ct2 = createCT2(Wf,fiberData);
                resolvedW = happinessCheck(ct2, [], ct2.coeffs, [], pref);
                if ~resolvedW
                    m(3) = 2*m(3)-1;
                end
            end
        end
        
        [~, tol] = getTol(Uf, pseudoLevel, tol, dom(2)-dom(1));
        [~, tol] = getTol(Vf, pseudoLevel, tol, dom(4)-dom(3));
        [~, tol] = getTol(Wf, pseudoLevel, tol, dom(6)-dom(5));
        
        %% Phase 3
        
        % Compute factor matrices
        [Q1,R1] = qr(Uf,0);
        I = DEIM(Q1);

        [Q2,R2] = qr(Vf,0);
        J = DEIM(Q2);

        [Q3,R3] = qr(Wf,0);
        K = DEIM(Q3);
        
        % with diagonal scaling
        D1 = diag(eps./max(min(abs(chebfun(Q1).coeffs)),eps));
        D2 = diag(eps./max(min(abs(chebfun(Q2).coeffs)),eps));
        D3 = diag(eps./max(min(abs(chebfun(Q3).coeffs)),eps));
        cf3f.cols = chebfun(Q1*D1, [dom(1),dom(2)], pref);
        cf3f.rows = chebfun(Q2*D2, [dom(3),dom(4)], pref);
        cf3f.tubes = chebfun(Q3*D3, [dom(5),dom(6)], pref);
        cf3f.core = invtprod(invtprod(evalTensor(I,J,K, ff,vectorize),Q1(I,:),Q2(J,:),Q3(K,:)),D1,D2,D3);
        
        % chebfun simplification
        cf3f.cols = simplify(cf3f.cols, pref.chebfuneps, 'globaltol');
        cf3f.rows = simplify(cf3f.rows, pref.chebfuneps, 'globaltol');
        cf3f.tubes = simplify(cf3f.tubes, pref.chebfuneps, 'globaltol');
        
    end
    
    % Sample Test
    if ( passSampleTest && restarts <= maxRestarts )
        cf3fhandle = @(x,y,z) cf3f.feval(x,y,z);
        happy = sampleTest(f, cf3fhandle, tol, dom);
    else
        happy = 1;
    end
    
    % Restart
    if ~happy
        if restarts + 1 == maxRestarts
            warning('chebfun3f: max number of restarts reached')
            return
        end
        
        % Increase n
        n(1) = floor(sqrt(2)^(floor(2*log2(n(1))) + 1)) + 1;
        n(2) = floor(sqrt(2)^(floor(2*log2(n(2))) + 1)) + 1;
        n(3) = floor(sqrt(2)^(floor(2*log2(n(3))) + 1)) + 1;
        restarts = restarts + 1;
        
        % Ensure r is large enougth for
        % (1,r,r) functions
        if r(1) > 1 || r(2) > 1 || r(2) > 1
            if r(1) < 3
                r(3) = max(6,2*r(3));
                r(2) = max(6,2*r(2));
            elseif r(2) < 2
                r(1) = max(6,2*r(1));
                r(3) = max(6,2*r(3));
            elseif r(3) < 2
                r(1) = max(6,2*r(1));
                r(2) = max(6,2*r(2));
            end
        end
        
        % very low-rank functions
        r(1) = max(r(1),3);
        r(2) = max(r(2),3);
        r(3) = max(r(3),3);
    end
end
end

%% TODO replace this function by the functionality in Chebfun
function F = Vals2ChebCoeffsMat(n)
% maps function evaluations at n Chebyshev nodes to Chebyshev coefficients
if n < 2
    warning('n too small')
end
F = zeros(n);
cheb = @(i,n) cos((i-1).*pi/(n-1));
xx = cheb(1:n,n);
T = @(i,x) cos((i-1)*acos(x));
for i = 1:n
    
    F(i,:) = T(i,xx);
    
end
F(:,1) = F(:,1)/2;
F(1,:) = F(1,:)/2;
F(:,n) = F(:,n)/2;
F(n,:) = F(n,:)/2;
F = (2/(n-1)).*F;
end

%% Additional Functions

%% Adaptive Cross Approximation with full pivoting
function [Ac, At, Ar, rowInd, colInd] = ACA(A, tol, maxIter)
Ac = [];
Ar = [];
At = [];
rowInd = [];
colInd = [];
Aoriginal = A;

for iter = 1:maxIter
    
    [error,I2] = max(abs(A(:)));
    if isempty(error) || error < tol
        Ac = Aoriginal(:,colInd);
        Ar = Aoriginal(rowInd,:)';
        At = Aoriginal(rowInd,colInd);
        return
    end
    
    [I,J] = ind2sub(size(A), I2);
    rowInd = [rowInd, I];
    colInd = [colInd, J];
    
    A = A-A(:,J)*A(I,:)./A(I,J);
end
Ac = Aoriginal(:,colInd);
Ar = Aoriginal(rowInd,:)';
At = Aoriginal(rowInd,colInd);
end

%% Discrete Empirical Interpolation 
function indices = DEIM(U)
indices = [];
[~, I] = max(abs(U(:,1)));
indices = [indices,I];
for l = 2:size(U,2)
    c = U(indices,1:(l-1)) \ U(indices,l);
    r = U(:,l) - U(:,1:(l-1))*c;
    [~, I] = max(abs(r));
    indices = [indices,I];
end
end

%% Create temporary chebtech2 object
function ct2 = createCT2(W,data)
data.vscale = max(abs(W(:)));
ct2 = chebtech2(W, data);
ct2.coeffs = sum(abs(ct2.coeffs), 2);
end

%% Get suitable tolerances as in chebfun3 (see https://github.com/chebfun/chebfun/issues/1491)
function [relTol, absTol] = getTol(M, pseudoLevel, tolOld,domDiff)
relTol = 2*size(M,1)^(4/5) * pseudoLevel;
vscale = max(abs(M(:)));
cheb = @(i,n) -cos((i-1).*pi/(n-1));
points = 1:size(M,1);
points = cheb(points, size(M,1));
gradNorms = zeros([1,size(M,1)]);
for i = 1:size(M,2)
    gradNorms(i) = max(abs(diff(M(:,i)) ./ diff(points)'));
end
gradNorms = max(gradNorms);
absTol = max(max(domDiff.*gradNorms), vscale) * relTol;
absTol = max([absTol, tolOld, pseudoLevel]);
end

%% Random initialization of indices by drawing one index in each subintervall of equal length
function X = initializeIndexRandomly(r, maxVal)
box = floor(maxVal/r);
X = [];
for i = 1:r
    val = i*box + randi(box,1);
    X = [X,val];
end
end

%% Evaluate X x_1 U x_2 V x_3 W
function X = tprod(X,U,V,W)
n = size(X);
m = [size(U,1),size(V,1),size(W,1)];
X = reshape(U*reshape(X,[n(1),n(2)*n(3)]),[m(1),n(2),n(3)]);
X = permute(reshape(V*reshape(permute(X,[2,1,3]),[n(2),m(1)*n(3)]),[m(2),m(1),n(3)]),[2,1,3]);
X = permute(reshape(W*reshape(permute(X,[3,2,1]),[n(3),m(2)*m(1)]),[m(3),m(2),m(1)]),[3,2,1]);
end

%% Evaluate X x_1 inv(U) x_2 inv(V) x_3 inv(W) using backslash
function X = invtprod(X,U,V,W)
n = size(X);
m = [size(U,1),size(V,1),size(W,1)];
X = reshape(U\reshape(X,[n(1),n(2)*n(3)]),[m(1),n(2),n(3)]);
X = permute(reshape(V\reshape(permute(X,[2,1,3]),[n(2),m(1)*n(3)]),[m(2),m(1),n(3)]),[2,1,3]);
X = permute(reshape(W\reshape(permute(X,[3,2,1]),[n(3),m(2)*m(1)]),[m(3),m(2),m(1)]),[3,2,1]);
end
