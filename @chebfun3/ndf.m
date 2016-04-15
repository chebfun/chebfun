function out = ndf(F)
%NDF   Number of degrees of freedom (parameters) needed to represent a 
% Chebfun3 object.

[rank_F1, rank_F2, rank_F3] = rank(F);

if ( isempty(F) )
    out = 0;
else
    s = rank_F1*rank_F2*rank_F3; % No of entries in the core tensor.
    for i = 1:rank_F1
        s = s + length(F.cols(:, i));
    end
    for i = 1:rank_F2
        s = s + length(F.rows(:, i));
    end
    for i = 1:rank_F3
        s = s + length(F.tubes(:, i));
    end
    out = s;
end

end