sim_title = "h-ffa-24";
outer_folder = "h-runs/new/";
data_folder = outer_folder + sim_title + "/DATA/";
plots_folder = data_folder + "plots/";

fname = data_folder + 'sup-input.json';
json_str = fileread(fname);
json_arr = jsondecode(json_str);
Constants;

pulse{1}.lambda0 = json_arr.laser.pulses(1).lambda0;
pulse{1}.omg0 = 2 * pi * const.cl / pulse{1}.lambda0;

try
    pulse{2}.lambda0 = json_arr.laser.pulses(2).lambda0;
    pulse{2}.omg0 = 2 * pi * const.cl / pulse{2}.lambda0;
catch 
    plot_spectral_energy =      false;
end

clearvars fname json_str;

% Path to SupData.m class
addpath("/gpfs/fs1/home/aellio18/Desktop/SUPER/super/analysis")

% Plotting defaults
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultfigurecolor',[1 1 1])
set(groot,'defaultAxesFontSize',18);

% -------------------------------------------------------------------------

Spectrum = SupData(data_folder + 'streak.h5');
Spectrum.axes(1).name = "$\omega$";
grid.omg = Spectrum.axes(1).data;
Spectrum.axes(1).data = ...
    Spectrum.axes(1).data / pulse{1}.omg0; % Normalize frequency
Spectrum.axes(1).units = "[$\omega_0$]";
Spectrum.axes(2).data = Spectrum.axes(2).data * 100;
Spectrum.axes(2).units = "[cm]";
grid.z = Spectrum.axes(2).data;

[~, i_zero] = min(abs(grid.omg/pulse{1}.omg0));

shifted_Spectrum = fftshift(Spectrum.data,1).';
max_spectrum = max(max(shifted_Spectrum));

% -------------------------------------------------------------------------

% Find the 2nd and 3rd harmonics
harm2 = 2 * pulse{1}.omg0;

try 
    harm3 = 3 * pulse{2}.omg0;
catch 
    harm3 = 3 * pulse{1}.omg0;
end

% Get the midpoint between them
harm_mid = (harm2 + harm3)/2;

% Find the nearest index
[~,omg1_index] = min(abs(pulse{1}.omg0 - grid.omg));
[~,omg2_index] = min(abs(pulse{2}.omg0 - grid.omg));
[~,harm_mid_index] = min(abs(harm_mid - grid.omg));

% Get index of the focal position, zf.
zf = json_arr.laser.pulses(1).zf * 100; % m -> cm
[~,zf_index] = min(abs(zf - Spectrum.axes(2).data));

% Calculate the SCG
% SCG_Spectrum = Calc_Supercontinuum_Spectra(shifted_Spectrum,zf_index);

u = size(shifted_Spectrum);
SCG_Spectrum = zeros(u);
init_spectra = shifted_Spectrum(1,:);
% figure; plot(init_spectra); title('Init Spectra');

for iz = 1:u(1)
    SCG_Spectrum(iz,:) = shifted_Spectrum(iz,:) - init_spectra;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate the symlog of the SCG Spectrum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_plot = Spectrum.axes(1).data(i_zero:harm_mid_index);
y_plot = Spectrum.axes(2).data;
z_plot = SCG_Spectrum(:, i_zero:harm_mid_index);

C = -30;
SCG_symlog = zeros(size(z_plot));
SCG_symlog = sign(z_plot).*(log10(1+abs(z_plot)/(10^C)));

max_SCG_spectrum = max(max(SCG_symlog));
orders = 1e-8;

figure; 
imagesc(x_plot, y_plot, SCG_symlog); 
axis xy; 
c = redblueTecplot();
colormap(c); %set according to your data
caxis([-max(abs(SCG_symlog(:))) max(abs(SCG_symlog(:)))]);
xlabel("$\omega [\omega_0]$"); 
ylabel("$z$ [cm]"); title('SCG Spectrum Symlog')
filename = plots_folder + sim_title + "_Spectrum_SCG_symlog_2.jpg";
colorbar;
saveas(gcf,filename);

% Save SCG_symlog plotting lineouts later
save(plots_folder + 'SCG_symlog','SCG_symlog');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create lineouts
%%%%%%%%
hffa24 = load('/gpfs/fs2/lle-scratch/aellio18/h-runs/new/h-ffa-24/DATA/plots/SCG_symlog.mat');
hffa24_new = load('/gpfs/fs2/lle-scratch/aellio18/h-runs/new/h-ffa-24/P_cr_3/DATA/plots/SCG_symlog.mat');
hffa26 = load('/gpfs/fs2/lle-scratch/aellio18/h-runs/new/h-ffa-26/DATA/plots/SCG_symlog.mat');

set(groot, 'defaultLineLineWidth', 1.5)
figure;
plot(x_plot,hffa24.SCG_symlog(8499,:)); hold on;
plot(x_plot,hffa24_new.SCG_symlog(8499,:)); hold on;
plot(x_plot,hffa26.SCG_symlog(8499,:)); hold off;
ylim([0,15])
legend("h-ffa-24, P/Pcr=0.037", "h-ffa-24, P/Pcr=3.7", "h-ffa-26, P/Pcr=3.7")
xlabel("$\omega [\omega_0]$"); 
title('SCG Spectrum Symlog')