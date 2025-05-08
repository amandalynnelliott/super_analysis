% F = scatteredInterpolant(theta_grid(:), omg_grid(:), Ef_trans(i_zero:end,:), 'linear', 'none');
% 
% [theta_target, omg_target] = meshgrid(theta_g, grid.omgpos);
% Eg = F(theta_target, omg_target);  % Eg: [Nomega x Ntheta]

theta = linspace(0,max_theta,grid.Nr);
[theta_A,omg_A] = meshgrid(theta,grid.omg);
kperp_A = omg_A/const.cl.*sin(theta_A);

Eg = interp2(kperp_grid,omg_grid,Ef_trans,kperp_A,omg_A,'linear',0);

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

% theta = linspace(0,max_theta,grid.Nr);
% [theta_A,omg_A] = meshgrid(theta,grid.omg);
% kperp_A = omg_A/const.cl.*sin(theta_A);
% 
% Etest_interp = interp2(kperp_grid,omg_grid,Etest,kperp_A,omg_A,'linear',0);
% Ef_interp = interp2(kperp_grid,omg_grid,Ef_trans,kperp_A,omg_A,'linear',0);
% 
% max_temp = max(real(Etest), [], 'all');
% min_temp = min(real(Etest), [], 'all');
% 
% figure; imagesc(grid.omg, grid.kperp, real(Etest)'); axis xy;
% set(gca, 'ColorScale','log'); set(gca, 'CLim', [min_temp, max_temp]); 
% xlabel("omega [rad/s]"); ylabel("k [1/m]"); title("Test Field - Before")
% figure; imagesc(grid.omg, theta, real(Etest_interp)'); axis xy;
% xlabel("omega [rad/s]");  ylabel("Angle [rad]"); title("Test Field - After")
% 
% max_tran = max(real(Ef_trans), [], 'all');
% min_tran = min(real(Ef_trans), [], 'all');
% 
% figure; imagesc(grid.THz, grid.kperp, real(Ef_trans)'); axis xy;
% set(gca, 'ColorScale','log'); set(gca, 'CLim', [min_tran, max_tran]); 
% xlabel("f [THz]"); ylabel("k [1/m]"); title("Final THz - Before")
% 
% max_tran = max(real(Ef_interp), [], 'all');
% min_tran = min(real(Ef_interp), [], 'all');
% 
% figure; imagesc(grid.THz, theta, real(Ef_interp)'); axis xy;
% set(gca, 'ColorScale','log'); set(gca, 'CLim', [min_tran, max_tran]); 
% xlabel("f [THz]");  ylabel("Angle [rad]"); title("Final THz - After")