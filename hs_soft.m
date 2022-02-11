function threshed = hs_soft(x,tau)

    threshed = max(abs(x)-tau,0);
    threshed = threshed.*sign(x);
return