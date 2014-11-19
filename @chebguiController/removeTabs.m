function str = removeTabs(str)
% Remove tabs from inputs
for k = 1:numel(str)
    idx = 1;
    strk = str{k};
    while ( ~isempty(idx) )
        idx = strfind(strk, double(9));
        strk(idx) = [];
    end
    str{k} = strk;
end
end