function y = tv3d(x, tau, alphay, alphaz)
    Demention = 3;
    gamma = 1;   
    thresh = sqrt(2) * 2 * Demention * tau * gamma;
    y = zeros(size(x), 'like', x);
    for axis = 1 : 3
        if axis == 3
            t_scale = alphaz;
        elseif axis ==1
            t_scale = alphay;
        else
            t_scale = 1;
        end
        y = y + iht3(ht3(x, axis, false, thresh*t_scale), axis, false);
        y = y + iht3(ht3(x, axis, true, thresh*t_scale), axis, true);
    end
    y = y / (2 * Demention);
return 





