%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the time derivative of current density, dJ/dt, for THz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Take product: n(t) * E(t), then take FT.
% filter_n_E = Filter to get THz frequencies
% dJ_dt = iFT{i*omega / (i*omega - nu_ei) * (e^2/m) * filter_n_E}

consts = const.echarge^2 / const.me;   % SI units
nu_en = 10 * 1e12 * 2 * pi;  % Electron-neutral collision freq, 10 THz

% -------------------------------------------------------------------------
% On-axis
% -------------------------------------------------------------------------

product_t = neOA.data .* EOA.data;  % Both arrays are (t,z)
product_f = fftshift(ifft(product_t,[],1),1);  % Transform to freq-domain.

% Filter data (neOA * EOA) to get only THz frequencies.
filtered_product_f = product_f .* full_filter';

% Variables to multiply to product before FFT.
freq_variables = 1i .* grid.omg / (1i .* grid.omg - nu_en) .* consts;
final_product_f = freq_variables .* filtered_product_f;

% Transform to time-domain.
dJ_dt = fft(final_product_f,[],1)';
clearvars constants nu_en product_t product_f filtered_product_f final_product_f;

figure;
imagesc(grid.t, grid.z, real(dJ_dt)); 
axis xy; 
c = redblueTecplot(); colormap(c); colorbar; 
% set(gca, 'ColorScale', 'log');
xlabel("t [s]"); ylabel("z [cm]"); 
title("dJ/dt on-axis")
filename = plots_folder + sim_title + "_EOA_dJ_dt.jpg";
saveas(gcf,filename);
clearvars filename dJ_dt;


% -------------------------------------------------------------------------
% Each z-step:
% -------------------------------------------------------------------------
% ---- neF and EF is (r,t,z)
if ~exist('neF','var')
    neF = SupData(data_folder + 'neF.h5');
end

mkdir(plots_folder + "dJdt/")
dJdt_full = zeros(size_EF);

% Calculate dJ_dt for each z step.
for iz = 1:size_EF(3) - 1
    time_EF = EF.data(:,:,iz);  % Both arrays are (r,t)
    time_neF = neF.data(:,:,iz);

    product_t = time_neF .* time_EF;
    product_f = fftshift(ifft(product_t,[],2),2);  % Transform to freq-domain.

    % Filter data (neF * EF) to get only THz frequencies.
    filtered_product_f = full_filter .* product_f;

    % Multiply by variables.
    final_product_f = freq_variables .* filtered_product_f;
    
    % Transform to time-domain.
    dJdt_full(:,:,iz) = fft(final_product_f,[],2);

end

clearvars time_EF time_neF filtered_product_f final_product_f freq_variables;

% Plot each z-step
if plot_arr.GIF
        max_dJdt = max(real(dJdt_full), [], 'all');
        min_dJdt = min(real(dJdt_full), [], 'all');
end
figure;
for iz = 1:size_EF(3) - 1
    imagesc(grid.t, grid.r, real(dJdt_full(:,:,iz))); 
    axis xy; 
    c = redblueTecplot(); colormap(c); colorbar; 
    if plot_arr.GIF
         set(gca, 'CLim', [-max_dJdt, max_dJdt]);
    end
    a = colorbar;
    a.Label.String = "dJ/dt";
    xlabel("t [s]"); ylabel("r [m]"); 
    title("z = " + num2str(EF.axes(3).data(iz) + EF.axes(3).units))
    filename = plots_folder + "dJdt/" + sim_title + "_" + ...
        num2str(iz) + "_dJdt.jpg";
    saveas(gcf,filename);
    clearvars filename;
end

clearvars max_dJdt min_dJdt dJdt_full iz;




