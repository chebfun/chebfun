function h = plus(f, g)
%+   Plus for CHEBFUN3 objects.

% Copyright 2017 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

if ( ~isa(f, 'chebfun3') ) % ??? + CHEBFUN3
    h = plus(g, f);
    
elseif ( isempty(g) || isempty(f) ) % CHEBFUN3 + []    
    % Return empty CHEBFUN3.
    h = chebfun3();
    
elseif ( isa(g, 'double') )           % CHEBFUN3 + DOUBLE
    % Convert g to a CHEBFUN3.
    g = chebfun3(g, f.domain);
    h = plus(f, g);
    
elseif ( ~isa(g, 'chebfun3') )          % CHEBFUN3 + ???
    error('CHEBFUN:CHEBFUN3:plus:unknown', ...
        ['Undefined function ''plus'' for input arguments of type %s ' ...
        'and %s.'], class(f), class(g));
    
else                                     % CHEBFUN3 + CHEBFUN3
    % Domain Check:
    if ( ~domainCheck(f, g) )
        error('CHEBFUN:CHEBFUN3:plus:domain', 'Inconsistent domains.');
    end
    
    % Check for zero CHEBFUN3 objects:
    if ( iszero(f) )
        h = g;
    elseif ( iszero(g) )
        h = f;
    else        
        % Add together two nonzero CHEBFUN3 objects:
        % Developer Note: If f+g has values close to zero, i.e., when 
        % computing f+g is subject to cancellation errors, the tolerance
        % computed in Chebfun3 constructor, used to stop adding new terms
        % to the BTD, will be roughly eps^2 which we cannot attain anyway.
        % See constructor/getTol3D. Here, the idea is to send an absolute 
        % tolerance to the constructor for computing f+g. In other words,
        % we do not expect accuracy better than the condition number of the
        % plus operation times CHEBFUN3EPS. Let's compute the condition 
        % number of the plus operation f + g, i.e., ||f|| + ||g|| / ||f+g||
        % We first compute a rough approximate vscale for f+g:
        m = 51; % size of sampling grid

        % trigtech and chebtech sample at different points. The else-part of the following 
        % statement ensures that f and g are sampled at the same points even if they have different techs.
        if (~isPeriodicTech(f) && ~isPeriodicTech(g) ) || (isPeriodicTech(f) && isPeriodicTech(g) )
            hVals = sample(f, m, m, m) + sample(g, m, m, m);
        else 
            [xxf, yyf, zzf] = chebpts3(m, m, m, f.domain);
            [xxg, yyg, zzg] = chebpts3(m, m, m, g.domain);
            hVals = feval(f, xxf, yyf, zzf) + feval(g, xxg, yyg, zzg);
        end
        hVscale = max(abs(hVals(:)));
        vscaleFG = [vscale(f), vscale(g)];
        kappa = sum(vscaleFG)/hVscale;
        % Developer Note: Since kappa is big, if g ~= -f, this helps us 
        % find out cases prone to cancellation errors, i.e., in this case 
        % we expect less from the constructor. 
        pref = chebfunpref(); prefStruct = pref.cheb3Prefs;
        eps = prefStruct.chebfun3eps;
        tol = eps*kappa;
        h = chebfun3(@(x, y, z) feval(f, x, y, z) + feval(g, x, y, z), ...
            f.domain, 'eps', tol);         
%         h = compressed_plus(f, g); % bypass the constructor and call the
%         compressed version
   end 
    
end

end

function h = compressed_plus(f, g)
% This is an alternative plus algorithm for CHEBFUN3 objects. It takes the
% structure of F and G into account and does NOT call the constructror.

% The algorithm is a generalization of CHEBFUN2 plus to 3D. If 
% A = ACore x_1 ACols x_2 ARows x_3 ATubes, and 
% B = BCore x_1 BCols x_2 BRows x_3 BTubes, then
% C := A + B has a _block diagonal_ core tensor formed by ACore and BCore, 
% and we have: CCols = [ACols, BCols], CRows = [ARows, BRows], and 
% CTubes = [ATubes, BTubes].
% This result in ranks being doubled and the following code tries to 
% compress the factors and core tensor.

% It's drawback is that it might be much slower than calling the constructor.

[fCore, fCols, fRows, fTubes] = tucker(f);
[gCore, gCols, gRows, gTubes] = tucker(g);
[r1f, r2f, r3f] = rank(f); 
[r1g, r2g, r3g] = rank(g); 

allCols = [fCols gCols]; 
[QCols, RCols] = qr(allCols);

allRows = [fRows gRows]; 
[QRows, RRows] = qr(allRows);

allTubes = [fTubes gTubes]; 
[QTubes, RTubes] = qr(allTubes);

allCores = zeros(r1f+r1g, r2f+r2g, r3f+r3g);
allCores(1:r1f, 1:r2f, 1:r3f) = fCore;
allCores(r1f+1:end, r2f+1:end, r3f+1:end) = gCore;
allCores = chebfun3.txm(chebfun3.txm(chebfun3.txm(allCores, RCols, 1), ...
    RRows, 2), RTubes, 3);
[core, U1, U2, U3] = chebfun3.discreteHOSVD(allCores);

% No truncation yet
% % Compress the format if possible. Warning: In contrast to 2D, truncating
% a 3D tensor (even in our Tucker format) does NOT give an optimal
% approximation, but it usually gives good approximations.
vf = vscale(f); 
vg = vscale(g);
vscl = max(vf, vg); 
% Truncate mode-1 singular values that fall below eps*vscale: 
idx1 = find( hosv(core,1) > eps * vscl, 1, 'last');
idx2 = find( hosv(core,2) > eps * vscl, 1, 'last');
idx3 = find( hosv(core,3) > eps * vscl, 1, 'last');
if ( isempty(idx1) || isempty(idx1) || isempty(idx1) )
    % Example: Lf - div(grad(f)) is zero for the function from Guide18.m
    h = chebfun3(0);
    return
end
QCols = QCols(:, 1:idx1);
U1 = U1(1:idx1, 1:idx1);

QRows = QRows(:, 1:idx2); 
U2 = U2(1:idx2, 1:idx2);

QTubes = QTubes(:, 1:idx3); 
U3 = U3(1:idx3, 1:idx3);
% or instead of truncating QCols and U1, first multiply them and truncate
% the resulting quasimatrix.

core = core(1:idx1, 1:idx2, 1:idx3);

% Form the output:
h = f;
h.cols = QCols*U1;
h.rows = QRows*U2;
h.tubes = QTubes*U3;
h.core = core;

end

%%
function sv = hosv(T,n)
% Mode-n singular values of an all-orthogonal tensor T.
sv = zeros(size(T,n),1);
if n == 1
    for i=1:size(T,n)
        sv(i) = norm(squeeze(T(i,:,:)), 'fro');
    end
elseif n == 2
    for i=1:size(T,n)
        sv(i) = norm(squeeze(T(:,i,:)), 'fro');
    end
elseif n == 3
    for i=1:size(T,n)
        sv(i) = norm(T(:,:,i), 'fro');
    end
end

end

%%
% f = chebfun(@(x) sin(x),[0,1]);
% g = chebfun(@(x) sqrt(1-cos(x).^2),[0,1]);
% h = f-g % fast, because coeffs are added not values. But here is what happens 
% % if we force 1D chebfun to call the constructor for doing plus:
% h = chebfun(@(x) f(x)-g(x),[0,1])
% % Warning: Function not resolved using 65537 pts. Have you tried 'splitting on'? 
% % > In chebfun/constructor>constructorNoSplit (line 120)
% %   In chebfun/constructor (line 63)
% %   In chebfun (line 219) 

% 2 related problems for Chebfun2:
% 1) computing laplacian of a harmonic function with "compression_plus" being
% disabled.
% 2) https://github.com/chebfun/chebfun/issues/1536
% Here we are calling 2D constructor for "times" on an object close to zero.
