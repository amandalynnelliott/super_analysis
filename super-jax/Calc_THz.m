%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the energy contained in the THz portion of the spectrum. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the energy in the zmin slice, over all frequencies. 
% Check that this energy is equal to the energy from the diagnostic,
% ---- U(zmin), U.data(0)
Uinit = U(1);
grid.THz = grid.omg / (2 * pi * 1e12);
U_THz_max = 30; % Maximum THz frequency to integrate to.

init_spectrum = streak(:,1);
init_spectral_energy = trapz(init_spectrum,1) * grid.domg;
energy_diff = abs(Uinit - init_spectral_energy);
% disp(['Energy difference: ' + energy_diff])

if energy_diff > 1e-3
    warning('THz energy calculation is incorrect.')
end

clearvars Uinit init_spectrum init_spectral_energy energy_diff;

% Now get energy, but integrated over a specific d_omega (only the THz
% range: [0, 100]THz)
% --- Do this for all z. 
% --- Plot THz energy vs z. 

[~,i_zero] = min(abs(grid.THz));
grid.THz_half = grid.THz(i_zero:end);
[~,i_temp] = min(abs(grid.THz_half - U_THz_max));

spectral_energy_THz = ...
    trapz(streak(1:i_temp,:),1) * grid.domg;

figure;
plot(grid.z,spectral_energy_THz);
xlabel("z [m]"); title('Spectral Energy of THz')
filename = plots_folder + "THz_spectral_energy.jpg";
saveas(gcf,filename);
clearvars filename spectral_energy_THz U_THz_max;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THz Field and dJ/dt Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate final THz field.
% -------------------------------------------------------------------------
% Constants for changing between spectral and real domains.
const.R2S = length(grid.t) * grid.dt;
const.S2R = grid.domg ./ (2*pi);

% Filter out non-THz frequencies.
freq_max = 30; 
freq_min = 0;
% freq_min = grid.domg / THz; % Make min freq one grid.omg from zero.

filter_gaussian = '1'; % 1: Hard, 2: Gaussian
filter_half = true;    % Half only cuts off at zero. False is true Gaussian.

full_filter = Calc_THz_filter(filter_gaussian, filter_half, freq_min, freq_max, grid.THz);
clearvars freq_max freq_min filter_gaussian filter_half;

% Efinal -> ffEfinal_spectral = fft(Efinal.data)
Efinal = EF.r(:,:,end) + 1j * EF.i(:,:,end);
Efinal_spectral = const.R2S .* ifftshift(ifft(Efinal,[],1),1);
Efinal_spectral_THz = full_filter' .* Efinal_spectral;
Efinal_THz_t = const.S2R .* fft(fftshift(Efinal_spectral_THz,1),[],1);

% Plot to make sure that first and second harmonics are there.
% figure; imagesc(grid.omg, grid.r,abs(Efinal_spectral)'); axis xy; 
% xlabel("$\omega$ [rad/s]"); ylabel("r [m]"); 
% title("Final E Field - Spectral - Before Filter");
% 
% figure; imagesc(grid.THz, grid.r,abs(Efinal_spectral_THz)'); axis xy;
% xlabel("f [THz]"); ylabel("r [m]");
% title("Final E Field - Spectral - After Filter");

max_THz = max(real(Efinal_THz_t), [], 'all');
min_THz = min(real(Efinal_THz_t), [], 'all');

figure; imagesc(grid.t, grid.r, real(Efinal_THz_t)'); 
axis xy; 
c = redblueTecplot();
colormap(c); 
colorbar; 
set(gca, 'CLim', [-max_THz, max_THz]);
a = colorbar;
a.Label.String = "THz Field";
xlabel("t [s]"); ylabel("r [m]"); title("Final THz Field")
filename = plots_folder + "Efinal_THz.jpg";
saveas(gcf,filename);
clearvars filename min_THz Efinal_spectral Efinal_spectral_THz;



% Calculate dJ/dt on axis.
% -------------------------------------------------------------------------

consts = const.echarge^2 / const.me;   % SI units
nu_en = 10 * 1e12 * 2 * pi;  % Electron-neutral collision freq, 10 THz

% -------------------------------------------------------------------------
% On-axis
% -------------------------------------------------------------------------

EOA_c = EOA.r + 1j * EOA.i;
% clearvars EOA;

product_t = neOA .* EOA_c;  % Both arrays are (t,z)
product_f = fftshift(ifft(product_t,[],1),1);  % Transform to freq-domain.

% Filter data (neOA * EOA) to get only THz frequencies.
filtered_product_f = product_f .* full_filter';

% Variables to multiply to product before FFT.
freq_variables = 1i .* grid.omg / (1i .* grid.omg - nu_en) .* consts;
final_product_f = freq_variables .* filtered_product_f;

% Transform to time-domain.
dJ_dt = fft(final_product_f,[],1);
% clearvars constants nu_en product_t product_f filtered_product_f final_product_f;

figure;
imagesc(grid.t, grid.z, real(dJ_dt)'); 
axis xy; 
c = redblueTecplot(); colormap(c); colorbar; 
% set(gca, 'ColorScale', 'log');
xlabel("t [s]"); ylabel("z [m]"); 
title("dJ/dt on-axis")
filename = plots_folder + "EOA_dJ_dt.jpg";
saveas(gcf,filename);
clearvars filename dJ_dt;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filter Spectrum to THz frequencies.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filtered_Spectrum = full_filter(i_zero:end)' .* streak;

figure; 
imagesc(grid.THz_half(1:i_temp),grid.z,...
    filtered_Spectrum(1:i_temp,:)'); 
axis xy; 
colorbar; 
% set(gca, 'ColorScale','log'); 
% set(gca, 'CLim', [orders*max_spectrum,max_spectrum]); % Orders of magnitude
colormap(parula);
xlabel("f [THz]"); 
ylabel("z [m]"); title('THz Spectrum')
filename = plots_folder + "Spectrum_THz.jpg";
saveas(gcf,filename);
clearvars filename filtered_Spectrum;