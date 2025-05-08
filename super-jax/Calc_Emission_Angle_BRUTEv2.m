% Smoothed version of Calc_Emission_Angle_BRUTE.m

% Initialize
Eg = zeros(size(theta_g,2), size(grid.omgpos,2));

% Set how tight to match theta
% (slightly larger delta_lim to capture multiple k_perp points)
delta_lim = 0.02; % You can tune this for smoothing effect

% Preallocate Eg
Eg = zeros(length(grid.omgpos), length(theta_g));

for iomg = 1:length(grid.omgpos)
    for itheta = 1:length(theta_g)

        % Find all k_perp indices whose theta is close enough
        delta = abs(theta_g(itheta) - theta_grid(iomg,:));
        idx_near = find(delta <= delta_lim);

        if ~isempty(idx_near)
            % Average contributions (could also use sum if you prefer)
            Eg(iomg, itheta) = mean(Ef_trans(iomg, idx_near));
        else
            Eg(iomg, itheta) = 0;
        end

    end
end

[~, ind] = min(abs(grid.THz - 30));
max_temp = max(abs(Eg(1:ind - i_zero,:)).^2, [], 'all');
min_temp = min(abs(Eg(1:ind - i_zero,:)).^2, [], 'all');

figure;
imagesc(grid.THz(i_zero:ind), theta_g, abs(Eg(1:ind - i_zero,:)').^2);
axis xy;
set(gca, 'ColorScale', 'linear');
set(gca, 'CLim', [min_temp, max_temp]); % Adjust depending on signal
colorbar;
xlabel("f [THz]");
ylabel("Angle [rad]");
title("Final THz - After");
% xlim([0,30]);


% Calculate the maximum emission angle (theta_C)
[max_val, max_idx] = max(real(Eg(:)));
[i_omega, i_theta] = ind2sub(size(Eg), max_idx);
theta_C = theta_g(i_theta);

fprintf('Peak emission at angle theta_C = %.4f radians\n', theta_C);