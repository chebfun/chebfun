function [M, S] = convertOperator(disc, item)
    % Developer note: In general we can't represent functional
    % blocks via coeffs. To get around this we instantiate a
    % TRIGCOLLOC discretization and convert it to coefficient space
    % using COEFFS2VALS(). (Note it's COEFFS2VALS() rather than
    % VALS2COEFFS() because it's a right-multiply (I think..).)

    % For convenience:
    dim = disc.dimension;
    dom = disc.domain;

    % Create a TECH and a VALSDISCRETIZATION:
    tech = disc.returnTech;
    tech = tech();
    valsDisc = tech.returnValsDisc;
    valsDisc = valsDisc(item, dim, dom);
    M = matrix(valsDisc);
    S = [];
    
    % Convert operator to and frmo coefficient space:
    M = tech.vals2coeffs(tech.coeffs2vals(M.').');

    % Try to make M sparse
    M_flat = M(:);
    tol = norm(M_flat, inf)*1e-10;        
    if ( norm(imag(M_flat), inf) < 10*tol )
        M = real(M);
    elseif ( norm(real(M_flat), inf) < 10*tol )
        M = imag(M);
    end      
    M(abs(M)<tol) = 0;
    if ( nnz(M)/numel(M) < .1 )
        M = sparse(M);
    end
    
end