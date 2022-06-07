function X = initializeIndexRandomly(r, maxVal)
box = floor(maxVal/r);
X = [];
for i = 1:r
    val = i*box + randi(box,1);
    X = [X,val];
end
end

