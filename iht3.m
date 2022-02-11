function y = iht3(w, ax, shift)
    s = size(w);
    y = zeros(s, 'like', w);
    C = 1 / sqrt(2);
    m = floor(s(ax) / 2);
    if ax == 1
        y(1:2:end, :, :) = C * (w(1:m, :, :) - w((m + 1):end, :, :));
        y(2:2:end, :, :) = C * (w(1:m, :, :) + w((m + 1):end, :, :));
    elseif ax == 2
        y(:, 1:2:end, :) = C * (w(:, 1:m, :) - w(:, (m + 1):end, :));
        y(:, 2:2:end, :) = C * (w(:, 1:m, :) + w(:, (m + 1):end, :));
    else
        y(:, :, 1:2:end) = C * (w(:, :, 1:m) - w(:, :, (m + 1):end));
        y(:, :, 2:2:end) = C * (w(:, :, 1:m) + w(:, :, (m + 1):end));
    end
    if shift
        y = circshift(y, 1, ax);
    end
return