function pass = domainCheck(f, g)

hs = max(hscale(f), hscale(g));
pass = norm(f.domain([1, end]) - g.domain([1, end]), inf) < 1e-15*hs;

end
