function out = normest(f)

out = 0;
for k = 1:numel(f.funs);
    out = out + normest(f.funs{k});
end

end