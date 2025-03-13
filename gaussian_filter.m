function [freq, filter] = gaussian_filter(f_range, f_center, bandwidth, n_points)
    
    % Validate inputs
    if numel(f_range) ~= 2 || f_range(1) >= f_range(2)
        error('f_range must be a 2-element vector with f_min < f_max.');
    end

    if f_center < f_range(1) || f_center > f_range(2)
        error('f_center must lie within f_range.');
    end

    f_min = f_range(1);
    f_max = f_range(2);
    freq = linspace(f_min, f_max, n_points);

    power = 10;

    filter = exp(-((freq - f_center).^power / (2 * bandwidth^power)));
end