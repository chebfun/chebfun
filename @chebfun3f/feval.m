function out = feval(cf3f,x,y,z)
    u = cf3f.cols(x);
    v = cf3f.rows(y);
    w = cf3f.tubes(z);
    
    r = cf3f.rank();
    
    out = cf3f.core;
    out = u*reshape(out,r(1),r(2)*r(3));
    out = v*reshape(out,r(2),r(3))*w';
end

