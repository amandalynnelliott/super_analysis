temp_size = size(theta_grid);
Eg = zeros(temp_size);

delta_lim = 0.01;

for iomg = 1:length(grid.omgpos)

    for itheta = 1:length(theta_g)

        delta = abs(theta_g(itheta) - theta_grid(iomg,:));  % Array of differences
        [min_delta, ik] = min(delta);  % Find the closest k_perp index

        if min_delta <= delta_lim
            Eg(iomg, itheta) = Ef_trans(iomg, ik);  % Only assign if close enough
        else
            Eg(iomg, itheta) = 0;  
        end

    end
end

[~, ind] = min(abs(grid.THz - 30));
max_temp = max(abs(Eg(1:ind - i_zero,:)).^2, [], 'all');
min_temp = min(abs(Eg(1:ind - i_zero,:)).^2, [], 'all');

theta_deg = theta_g * 180/pi;

figure;
surf(grid.THz(i_zero:ind), theta_deg, abs(Eg(1:ind - i_zero + 1,:)').^2, 'EdgeColor', 'none');
view(2);
set(gca, 'ColorScale', 'linear');
set(gca, 'CLim', [min_temp, max_temp]); % Adjust depending on signal
colorbar;
xlabel("f [THz]");
ylabel("Angle [degrees]");
title("Final THz - After");
% xlim([0,30]);

