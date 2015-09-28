function h = plus(f, g)
%+   Plus for SPHEREFUN objects.
%
% F + G adds F and G. F and G can be scalars or SPHEREFUN objects.

if ( ~isa(f, 'spherefun') ) % ??? + SPHEREFUN
    
    h = plus(g, f);
    
elseif ( isempty(g) ) % SPHEREFUN + []
    
    h = f; 
    
elseif ( isempty(f) ) % [] + SPHEREFUN
    
    h = g; 
    
elseif ( isa( g, 'double' ) )           % SPHEREFUN + DOUBLE
    
    g = compose( 0*f,@plus, g);   % promote double to object class.  
    h = plus(f, g); 
    
elseif ( ~isa(g, 'spherefun') )          % SPHEREFUN + ???
    
    error( 'SPHEREFUN:plus:unknown', ...
        ['Undefined function ''plus'' for input arguments of type %s ' ...
        'and %s.'], class( f ), class( g ));
        
else                                     % SPHEREFUN + SPHEREFUN
    
    % Domain Check:
    if ( ~domainCheck(f, g) )
        error('SPHEREFUN:plus:domain', 'Inconsistent domains.');
    end
    
    % Check for zero SPHEREFUN objects:
    if ( iszero(f) )
        h = g;
    elseif ( iszero(g) )
        h = f;
    else
        % Add together two nonzero SPHEREFUN objects:
        % The algorithm is as follows: Split f and g into their plus/minus
        % components.  Do the compression_plus algorithm described in
        % @separableApprox/compression plus on each pair of plus and minus
        % components.
        
        % Check if the pole is non-zero.  If it is then we need to strip
        % out the column and row that deal with this in the plus piece
        % before doing compression plus.  The reason is that what we feed
        % compression plus needs to be zero at the poles if what is to be
        % returned is zero at the poles.  Including one column that is
        % non-zero at the poles can screw everything up.
        [f,fPole] = extractPole(f);
        [g,gPole] = extractPole(g);

        [fp,fm] = partition(f);
        [gp,gm] = partition(g);
        
        hp = plus@separableApprox(fp,gp);
        r = size(hp.cols,2);
        hp.idxPlus = 1:r;
        % Indices or locations of the pivots do not make sense after a
        % compression plus, so we set them to NaN.
        hp.pivotLocations = nan(r,2);  % This should be done at the separableApprox level.
        
        hm = plus@separableApprox(fm,gm);
        r = size(hm.cols,2);
        hm.idxMinus = 1:r;
        hm.pivotLocations = nan(r,2);
                
        if ( ~isempty( fPole ) ) || ( ~isempty( gPole ) )
            % Set tolerance for determining if fPole+gPole=0.
            tol = eps*max(vscale(f),vscale(g)); 
            g = addPoles(fPole,gPole,tol);
            % Handle the rare case that g is zero and hp is not empty
            if ( g.pivotValues == 0 ) && ( ~isempty( hp ) )
                % Set g to empty spherefun
                g = spherefun([]);
            end
            hp.cols = [g.cols hp.cols];
            hp.rows = [g.rows hp.rows];
            hp.pivotValues = [g.pivotValues;hp.pivotValues];
            hp.pivotLocations = [g.pivotLocations;hp.pivotLocations];
            hp.idxPlus = 1:size(hp.cols,2);
            hp.nonZeroPoles = ~isempty(g);
        end
            
        % Put pieces back together.
        h = combine(hp,hm);

    end 
    
end

end

function f = addPoles(f,g,tol)

if isempty(g)
    return;
elseif isempty(f)
    f = g;
    return;
end

fmean = mean(f.rows);
gmean = mean(g.rows);

cols = (fmean/f.pivotValues)*f.cols + (gmean/g.pivotValues)*g.cols;

% If cols is numerically zero then return an empty result
if norm(cols) <= tol
    fmean = 0;
    gmean = 0;
    pivot = 0;
    nonZeroPoles = 0;
else
    pivot = 1;
    nonZeroPoles = 1;
end

f.cols = (fmean/f.pivotValues)*f.cols + (gmean/g.pivotValues)*g.cols;
f.rows = chebfun('1',f.domain(1:2),'trig');
f.pivotValues = pivot;
% No idea what indices or locations should be after plus
f.pivotLocations = [nan nan];
f.nonZeroPoles = nonZeroPoles;

end



