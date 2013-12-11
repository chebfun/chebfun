function out = numColumns(f)

if ( isempty(f) )
    out = 0;
elseif ( numel(f) == 1 )
    out = size(f.funs{1}, 2);
else
    out = numel(f);
end

end