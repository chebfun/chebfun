function scl = scale( disc, uFun ) 

for j = 1:numel( uFun ) 
    tmp = chebtech2.coeffs2vals( uFun{j} );
    scl(j) = max( abs(tmp) );
end

scl = max(scl);
end