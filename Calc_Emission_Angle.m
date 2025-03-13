%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Keys
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sim_title = "h-ga-26";
outer_folder = "h-runs/new/";
data_folder = outer_folder + sim_title + "/DATA/";
plots_folder = data_folder + "plots/";

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
%% Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load electric field: E(r,t,z)
Constants;
EF = SupData(data_folder + 'EF-re.h5',true);
grid.r = EF.axes(1).data;
grid.t = EF.axes(2).data;
grid.z = EF.axes(3).data;

% Load input array
fname = data_folder + 'sup-input.json';
json_str = fileread(fname);
json_arr = jsondecode(json_str);

clearvars json_str fname;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Quantities for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pulse{1}.lambda0 = json_arr.laser.pulses(1).lambda0;
pulse{1}.omg0 = 2 * pi * const.cl / pulse{1}.lambda0;
try
    pulse{2}.lambda0 = json_arr.laser.pulses(2).lambda0;
    pulse{2}.omg0 = 2 * pi * const.cl / pulse{2}.lambda0;
catch 
    plot_spectral_energy =      false;
end

Spectrum = SupData(data_folder + 'streak.h5');
Spectrum.axes(1).name = "$\omega$";
grid.omg_rad = Spectrum.axes(1).data;
Spectrum.axes(1).data = Spectrum.axes(1).data / pulse{1}.omg0; % Normalize frequency
Spectrum.axes(1).units = "[$\omega_0$]";
grid.omg = Spectrum.axes(1).data;

[~, i_zero] = min(abs(grid.omg));

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
[~,omg1_index] = min(abs(pulse{1}.omg0 - grid.omg_rad));
[~,omg2_index] = min(abs(pulse{2}.omg0 - grid.omg_rad));
[~,harm_mid_index] = min(abs(harm_mid - grid.omg_rad));

% Get index of the focal position, zf.
zf = json_arr.laser.pulses(1).zf;
[~,zf_index] = min(abs(zf - Spectrum.axes(2).data));
zf_index = ceil(zf_index / json_arr.diag.fullskip);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prelim quantities for transformations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ell = 0; % orbital angular momentum of pulse --> Bessel of order 0
grid.rmax = json_arr.grid.xmax(2); % Get maximum R
grid.nr = json_arr.grid.nx(2); % Get number of grid points in r, Nr

grid.nt = json_arr.grid.nt;
grid.dt = EF.axes(2).dax;

[grid.r, grid.kperp] = DHTarrays(grid.nr,grid.rmax,ell);
[grid.T_R2K,grid.T_K2R] = DHTtransformMatrix(grid.nr,grid.rmax,ell);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transformations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Perform FFT for E(r,t,z) --> E(r,omg,z) 
% Perform Hankel transform for E(r,omg,z) --> E(k_perp,omg,z)
% -------------------------------------------------------------------------
EAp = grid.nt * grid.dt * ifftshift(ifft(EF.data(:,:,1),[],2),2);
EAp = grid.T_R2K*EAp;

EAp_zf = grid.nt * grid.dt * ifftshift(ifft(EF.data(:,:,zf_index),[],2),2);
EAp_zf = grid.T_R2K*EAp_zf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting transformed E-field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orders = 1e-6;

max_EAp = max(max(abs(EAp)));
% Plot E(k_perp, omg, z=0) colormap on k_perp vs omg 
figure; imagesc(grid.omg(i_zero:harm_mid_index), grid.kperp, abs(EAp(:,i_zero:harm_mid_index))); 
axis xy; 
colorbar; colormap(parula);
colorbar; set(gca, 'ColorScale','log'); 
set(gca, 'CLim', [orders*max_EAp,max_EAp]); % Orders of magnitude
ylabel("$k_\perp [m^{-1}]$"); xlabel("$\omega [\omega_0]$");  title("$E(z=0)$")
filename = plots_folder + sim_title + "_EF-transformed_z0.jpg";
saveas(gcf,filename);
clearvars filename;

max_EAp = max(max(abs(EAp_zf)));
% Plot E(k_perp, omg, z=0) colormap on k_perp vs omg 
figure; imagesc(grid.omg(i_zero:harm_mid_index), grid.kperp, abs(EAp_zf(:,i_zero:harm_mid_index))); 
axis xy; 
colorbar; colormap(parula);
colorbar; set(gca, 'ColorScale','log'); 
set(gca, 'CLim', [orders*max_EAp,max_EAp]); % Orders of magnitude
ylabel("$k_\perp [m^{-1}]$"); xlabel("$\omega [\omega_0]$");  title("$E(z=z_f)$")
filename = plots_folder + sim_title + "_EF-transformed_zf.jpg";
saveas(gcf,filename);
clearvars filename;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate angle of emmission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% At what angle are the sidebands emitted? 
% By changing vf of ff pulses, do we cange theta of SCG?


% Calculate k_z = sqrt(omg^2 * n(omg) /c^2 -kperp^2) 
% grid.kz = sqrt(grid.omg2D.^2)

% EXTRAPOLATE FOR OTHER GAS ATOMS
%Real refractive index at 1 atm and 273 K, 1.785e-4 gm/cm^3 (can be
%scaled by mass density) from: https://refractiveindex.info/?shelf=main&book=He&page=Ermolov
scale = 1; % Put density fraction here!
lamb = (2.0*pi*const.cl./grid.omg_rad)/1e-6;

dindex = scale*( 2.50141e-3./(91.012 - lamb.^(-2)) ...
                + 5.00283e-4./(87.892 - lamb.^(-2)) ...
                + 5.22343e-2./(214.02 - lamb.^(-2)) ); 

nindex = 1 + dindex;

temp_omg = grid.omg_rad .* nindex / const.cl;

% Define meshgrid for k_perp and omg
[grid.kperp2D,grid.omg2D] = meshgrid(grid.kperp, temp_omg);
grid.kz = sqrt((grid.omg2D.^2) - grid.kperp2D.^2);

% Calculate theta = atan(k_perp/k_z)
theta = atan(grid.kperp2D./grid.kz);

% Plot theta vs omg
figure; imagesc(abs(theta)); 
axis xy; 
colorbar; colormap(parula);
% colorbar; set(gca, 'ColorScale','log'); 
ylabel("$N_r$"); xlabel("$N_t$");  
% ylabel("$r [m]$"); xlabel("$t [s]$");  
title("$\theta$")
filename = plots_folder + sim_title + "_theta.jpg";
saveas(gcf,filename);
clearvars filename;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r,k] = DHTarrays(N,R,ell)

    %This function defines the r and k array for a Hankel Transform
    %N is the total number of grid points
    %R is the maximum radius
    %ell is the order of teh Bessel function
    
    filename = ['BesselZeroes',num2str(ell),'.csv'];
    Zeroes = load(filename);
    r = R*Zeroes(1:N)/Zeroes(N);
    k = Zeroes(1:N)/R;

end

function [T_R2K,T_K2R] = DHTtransformMatrix(N,R,ell)

    %This function defines the Transform Matrix for the Hankel Transform
    %N is the total number of grid points
    %R is the maximum radius
    %ell is the order of the Bessel function
    %T_R2K transforms from r space to k space
    %T_K2R transforms from k space to r space
    
    filename = ['BesselZeroes',num2str(ell),'.csv'];
    Zeroes = load(filename);
    alpha = Zeroes(1:N);
    
    kR = alpha(N);
    kmax = kR/R;
    T = besselj(ell,alpha*alpha'/kR) ./ ( abs(besselj(ell+1,alpha))*abs(besselj(ell+1,alpha))' ) / kR;
    
    %A = repmat(abs(besselj(1,alpha)),[1,N]);
    A = abs(besselj(ell+1,alpha));
    
    T_R2K = 4.0*pi*(R/kmax)*A.*T./A.';
    T_K2R = (1/pi)*(kmax/R)*A.*T./A.'; 

end