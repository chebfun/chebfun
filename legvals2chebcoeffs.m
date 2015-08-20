function c_cheb = legvals2chebcoeffs(v_leg)

c_leg = chebfun.idlt(v_leg);
c_cheb = leg2cheb(c_leg);

end