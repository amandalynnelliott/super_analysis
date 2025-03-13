%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the energy contained in the THz portion of the spectrum. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate the energy in the zmin slice, over all frequencies. 
% Check that this energy is equal to the energy from the diagnostic,
% ---- U(zmin), U.data(0)
Uinit = U.data(1);
grid.THz = grid.omg / (2 * pi * 1e12);
U_THz_max = 100; % Maximum THz frequency to integrate to.

init_spectrum = shifted_Spectrum(1,:);
init_spectral_energy = trapz(init_spectrum,2) * grid.domg;
energy_diff = abs(Uinit - init_spectral_energy);
disp(['Energy difference: ' + energy_diff])

if energy_diff > 1e-3
    warning('THz energy calculation is incorrect.')
end

clearvars Uinit init_spectrum init_spectral_energy energy_diff;

% Now get energy, but integrated over a specific d_omega (only the THz
% range: [0, 100]THz)
% --- Do this for all z. 
% --- Plot THz energy vs z. 

[~,i_zero] = min(abs(grid.THz));
[~,i_temp] = min(abs(grid.THz - U_THz_max));

spectral_energy_THz = ...
    trapz(shifted_Spectrum(:, i_zero:i_temp),2) * grid.domg;

figure;
plot(grid.z,spectral_energy_THz);
xlabel("$z$ [cm]"); title('Spectral Energy of THz')
filename = plots_folder + sim_title + "_THz_spectral_energy.jpg";
saveas(gcf,filename);
clearvars filename spectral_energy_THz U_THz_max;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get the THz temporal profile from the E-field 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants for changing between spectral and real domains.
const.R2S = length(grid.t) * grid.dt;
const.S2R = grid.domg ./ (2*pi);

% Filter out non-THz frequencies.
freq_max = 150; 
freq_min = 0;
% freq_min = grid.domg / THz; % Make min freq one grid.omg from zero.

filter_gaussian = '2'; % 1: Hard, 2: Gaussian
filter_half = true;    % Half only cuts off at zero. False is true Gaussian.

full_filter = Calc_THz_filter(filter_gaussian, filter_half, freq_min, freq_max, grid.THz);
clearvars freq_max freq_min filter_gaussian filter_half;

% Efinal -> ffEfinal_spectral = fft(Efinal.data);
Efinal_spectral = const.R2S .* fftshift(ifft(Efinal.data,[],2),2);

% Plot to make sure that first and second harmonics are there.
% figure; imagesc(grid.omg, grid.r,abs(Efinal_spectral)); axis xy; 
% xlabel("$\omega$ [rad/s]"); ylabel("r [m]"); 
% title("Final E Field - Spectral - Before Filter");

Efinal_spectral_THz = Efinal_spectral .* full_filter;

% Plot after filtering
% figure; imagesc(grid.omg,grid.r,abs(Efinal_spectral_THz)); axis xy;
% xlabel("$\omega$ [rad/s]"); ylabel("r [m]"); 
% title("Final E Field - Spectral - After Filter");

Efinal_THz_t = const.S2R .* fft(Efinal_spectral_THz,[],2);

max_THz = max(real(Efinal_THz_t), [], 'all');
min_THz = min(real(Efinal_THz_t), [], 'all');

figure; imagesc(grid.t, grid.r, real(Efinal_THz_t)); 
axis xy; 
c = redblueTecplot();
colormap(c); 
colorbar; 
set(gca, 'CLim', [-max_THz, max_THz]);
a = colorbar;
a.Label.String = "THz Field";
xlabel("t [s]"); ylabel("r [m]"); title("Final THz Field")
filename = plots_folder + sim_title + "_Efinal_THz.jpg";
saveas(gcf,filename);
clearvars filename max_THz min_THz Efinal_spectral Efinal_spectral_THz Efinal_THz_t;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filter Spectrum to THz frequencies.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filtered_Spectrum = shifted_Spectrum .* full_filter;

figure; 
imagesc(grid.THz(i_zero:i_temp),grid.z,...
    filtered_Spectrum(:, i_zero:i_temp)); 
axis xy; 
colorbar; 
% set(gca, 'ColorScale','log'); 
% set(gca, 'CLim', [orders*max_spectrum,max_spectrum]); % Orders of magnitude
colormap(parula);
xlabel("f [THz]"); 
ylabel("$z$ [cm]"); title('THz Spectrum')
filename = plots_folder + sim_title + "_Spectrum_THz.jpg";
saveas(gcf,filename);
clearvars filename filtered_Spectrum;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get THz at each z-step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate and plot the THz E-field in time at various positions in z.
if plot_arr.THz_full

    if exist('EF', 'var') ~= 1
        EF = SupData(data_folder + 'EF-re.h5',true);
    end
   
    mkdir(plots_folder + "THz_real/")
    clearvars size_EF iz;

    size_EF = size(EF.data);
    THz_full = zeros(size(EF.data));
    
    % Calculate the THz for each z.
    for iz = 1:size_EF(3) - 1
        temp_time = EF.data(:,:,iz);
        temp_freq = const.R2S .* fftshift(ifft(temp_time,[],2),2);
   
        temp_freq_filtered = temp_freq .* full_filter;
        THz_full(:,:,iz) = const.S2R .* fft(temp_freq_filtered,[],2);
    end
    clearvars iz temp_time temp_freq temp_freq_filtered;
        
    % Plotting
    if plot_arr.GIF
        max_THz = max(real(THz_full), [], 'all');
        min_THz = min(real(THz_full), [], 'all');
    end
    figure;
    for iz = 1:size_EF(3) - 1
        imagesc(grid.t, grid.r, real(THz_full(:,:,iz))); 
        axis xy; 
        c = redblueTecplot();
        colormap(c); 
        colorbar; 
        if plot_arr.GIF
             set(gca, 'CLim', [-max_THz, max_THz]);
        end
        a = colorbar;
        a.Label.String = "THz Field";
        xlabel("t [s]"); ylabel("r [m]"); 
        title("z = " + num2str(EF.axes(3).data(iz) + EF.axes(3).units))
        filename = plots_folder + "THz_real/" + sim_title + "_" + ...
            num2str(iz) + "_THz_real.jpg";
        saveas(gcf,filename);
        clearvars filename;
    end

    clearvars iz max_THz min_THz iz THz_full;
end
