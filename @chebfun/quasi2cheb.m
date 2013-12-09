function F = quasi2cheb(F)

if ( numel(F) < 2 )
    return
end

F = restrict(F, get(F, 'domain'));
F = num2cell(F);
F = [F{:}];

end