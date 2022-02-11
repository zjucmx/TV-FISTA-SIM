function w = ht3(x, ax, shift, thresh)
    s = size(x);
    w = zeros(s, 'like', x);
    C = 1 / sqrt(2);
    if shift
        x = circshift(x, -1, ax);
    end
    m = floor(s(ax) / 2);
    if ax == 1
        w(1:m, :, :) = C * (x(2:2:end, :, :) + x(1:2:end, :, :));
        w((m + 1):end, :, :) = hs_soft(C * (x(2:2:end, :, :) - x(1:2:end, :, :)), thresh);
        %w((m + 1):end, :, :) = hs_soft(w((m + 1):end, :, :), thresh);
    elseif ax == 2
        w(:, 1:m, :) = C * (x(:, 2:2:end, :) + x(:, 1:2:end, :));
        w(:, (m + 1):end, :) = C * (x(:, 2:2:end, :) - x(:, 1:2:end, :));
        w(:, (m + 1):end, :) = hs_soft(w(:, (m + 1):end, :), thresh);
    else
        w(:, :, 1:m) = C * (x(:, :, 2:2:end) + x(:, :, 1:2:end));
        w(:, :, (m + 1):end) = C * (x(:, :, 2:2:end) - x(:, :, 1:2:end));
        w(:, :, (m + 1):end) = hs_soft(w(:, :, (m + 1):end), thresh);
    end
return

