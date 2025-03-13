%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Keys
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sim_title = "THz-f50-20-5d";
outer_folder = "Jeremy/two_freq/";
data_folder = outer_folder + sim_title + "/DATA/";
plots_folder = data_folder + "plots/";
mkdir(plots_folder)

% To Plot?
plot_THz =                  true;    % Run THz post-processing diagnostic. 
plot_THz_full =             true;    % Plot THz filtered E(z) at each z location

% Scaling for spectra
orders =                    1e-5;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Constants;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load from Input File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fname = data_folder + 'sup-input.json';
json_str = fileread(fname);
json_arr = jsondecode(json_str);


pulse{1}.lambda0 = json_arr.laser.pulses(1).lambda0;
pulse{1}.omg0 = 2 * pi * const.cl / pulse{1}.lambda0;

try
    pulse{2}.lambda0 = json_arr.laser.pulses(2).lambda0;
    pulse{2}.omg0 = 2 * pi * const.cl / pulse{2}.lambda0;
catch 
    plot_spectral_energy =      false;
end

clearvars fname json_str;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load plotting class and preferences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Path to SupData.m class
addpath("/gpfs/fs1/home/aellio18/Desktop/SUPER/super/analysis")

% Plotting defaults
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultfigurecolor',[1 1 1])
set(groot,'defaultAxesFontSize',18);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Electric Field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % Initial E-field ---------------------------------------------------------
% Einit = SupData(data_folder + 'Einit-re.h5',true); 
% 
% figure; imagesc(Einit.axes(2).cax, Einit.axes(1).cax, abs(Einit.data)); 
% axis xy; 
% colorbar; colormap(parula);
% xlabel("t [s]"); ylabel("r [m]"); title("Initial Field")
% filename = plots_folder + sim_title + "_Einit.jpg";
% saveas(gcf,filename);
% clearvars filename;
% 
% 
% % Final E-field -----------------------------------------------------------
% Efinal = SupData(data_folder + 'Efinal-re.h5',true); 
% 
% figure; imagesc(Efinal.axes(2).cax, Efinal.axes(1).cax, abs(Efinal.data)); 
% axis xy; 
% colorbar; colormap(parula);
% xlabel("t [s]"); ylabel("r [m]"); title("Final Field")
% filename = plots_folder + sim_title + "_Efinal.jpg";
% saveas(gcf,filename);
% clearvars filename;
% 
% % E-field at each z step --------------------------------------------------
% EF = SupData(data_folder + 'EF-re.h5',true);


% Get arrays
Spectrum = SupData(data_folder + 'streak.h5');
Spectrum.axes(1).name = "$\omega$";
grid.omg = Spectrum.axes(1).data;
grid.r = EF.axes(1).data;
grid.dr = EF.axes(1).dax;
grid.t = EF.axes(2).data;
grid.dt = EF.axes(2).dax;
grid.z = EF.axes(3).data;
grid.dz = EF.axes(3).dax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_THz
    Calc_THz;
end

