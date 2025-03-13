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

figure; 
subplot(2,2,1)
plot(grid.omg, init_spectra); hold on; 
plot(grid.omg, shifted_Spectrum(zf_index,:)); hold on;
plot(grid.omg, SCG_Spectrum(zf_index,:)); hold off;
set(gca, 'YScale','log')
legend('z=0, init Spectra', 'z=zf, before', 'z=zf, after');
title('Log Scale');

subplot(2,2,2)
plot(grid.omg, SCG_Spectrum(zf_index,:)); 
set(gca, 'YScale','linear')
legend('z=zf, after');
title('Linear');

subplot(2,2,3)
plot(grid.omg, SCG_Spectrum(zf_index,:));
symlog('y',-32)
% set(gca, 'YScale','log')
legend('z=zf, after');
title('Symlog');

% subplot(2,2,3)
% plot(grid.omg, init_spectra); hold on;
% plot(grid.omg, shifted_Spectrum(zf_index,:)); hold on;
% plot(grid.omg, SCG_Spectrum(zf_index,:)); hold off;
% set(gca, 'YScale','linear')
% legend('z=0, init Spectra', 'z=zf, before', 'z=zf, after');
% title('Linear');

subplot(2,2,4)
plot(grid.z, SCG_Spectrum(:,omg1_index)); hold on;
plot(grid.z, SCG_Spectrum(:,omg2_index)); hold off;
legend('omg1', 'omg2');
title('Linear');


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
filename = plots_folder + sim_title + "_Spectrum_SCG_symlog.jpg";
saveas(gcf,filename);

% x = linspace(-50,50,1e4+1);
% y1 = x;
% y2 = sin(x);
% y3 = x - sin(x);
% 
% plot(x,y1,x,y2,x,y3);
% symlog(gca, 'xy', -1.7);
% 
% Example:
% x = linspace(-50,50,1e4+1);
% y1 = x;
% y2 = sin(x);
% 
% subplot(2,4,1)
% plot(x,y1,x,y2)
% 
% subplot(2,4,2)
% plot(x,y1,x,y2)
% set(gca,'XScale','log') % throws warning
% 
% subplot(2,4,3)
% plot(x,y1,x,y2)
% set(gca,'YScale','log') % throws warning
% 
% subplot(2,4,4)
% plot(x,y1,x,y2)
% set(gca,'XScale','log','YScale','log') % throws warning
% 
% subplot(2,4,6)
% plot(x,y1,x,y2)
% symlog('x')
% 
% s = subplot(2,4,7);
% plot(x,y1,x,y2)
% symlog(s,'y') % can but don't have to provide s.
% 
% subplot(2,4,8)
% plot(x,y1,x,y2)
% symlog() % no harm in letting symlog operate in z axis, too.