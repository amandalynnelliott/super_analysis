%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Keys
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sim_title = "THz-f50-20-4b";
% outer_folder = "Jeremy/two_freq/";
% data_folder = outer_folder + sim_title + "/DATA/";
% plots_folder = data_folder + "plots/";

outer_folder = "../super-jax/output/";
sim_title = "Tanner-f10_ff_2color_extended_fixedU/";
% sim_title = "Tanner-f85-01b_gaus_streak2/";
% sim_title = "Tanner-f10_ff_2color_phase/";
data_folder = outer_folder + sim_title;
plots_folder = data_folder + "plots/";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load plotting class and preferences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path to SupData.m class
% addpath("/gpfs/fs1/home/aellio18/Desktop/SUPER/super/analysis")
addpath("/global/u1/a/aellio/super/analysis")

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

% For focusing on supercontinuum
grid.omg_short = grid.omg_rad(i_zero+20: harm_mid_index);

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
%% Calculate angle of emmission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRAPOLATE FOR OTHER GAS ATOMS -- only Argon so far
%Real refractive index at 1 atm and 273 K, 1.785e-4 gm/cm^3 (can be
%scaled by mass density) from: https://refractiveindex.info/?shelf=main&book=He&page=Ermolov

medium.atomsPmol = 1;

scale = json_arr.medium.atomdensity / const.ng0 / medium.atomsPmol;

lamb = (2.0*pi*const.cl./grid.omg_short)/1e-6;

dindex = scale*( 2.50141e-3./(91.012 - lamb.^(-2)) ...
                + 5.00283e-4./(87.892 - lamb.^(-2)) ...
                + 5.22343e-2./(214.02 - lamb.^(-2)) ); 

nindex = 1 + dindex;
%--------------------------------------------------------------------------

EAp_interp = abs(EAp(:,i_zero+20:harm_mid_index))';
EAp_zf_interp = abs(EAp_zf(:,i_zero+20:harm_mid_index))';

omg_dispersion = grid.omg_short .* nindex;

% Create (omega, theta) meshgrid.
[kperp_grid, omg_grid] = meshgrid(grid.kperp, grid.omg_short);
%[~,omg_grid_disp] = meshgrid(grid.kperp, omg_dispersion);
%theta_grid = asin(const.cl * kperp_grid ./ omg_grid_disp);
theta_grid = asin(const.cl * kperp_grid ./ omg_grid);

Etest = exp(-(theta_grid - 0.05).^2/(.02)^2);

theta = linspace(0,.1,100);
[theta_A,omg_A] = meshgrid(theta,grid.omg_short);
kperp_A = omg_A/const.cl.*sin(theta_A);

%EAp_zf_new = interp2(kperp_grid,omg_grid,EAp_zf_interp,kperp_A,omg_A,'linear',0);
EAp_zf_new = interp2(kperp_grid,omg_grid,Etest,kperp_A,omg_A,'linear',0);
%EAp_zf_new(isnan(EAp_zf_new)) = 0;
figure; imagesc(EAp_zf_new); axis xy;

a=1;

% valid_idx = (const.cl * kperp_grid ./ omg_grid_disp) <= 1;
% theta_grid(~valid_idx) = NaN;

theta_new = linspace(max(min(theta_grid(:)),-1), min(max(theta_grid(:)),1), 512);
[theta_new_grid, omg_new_grid] = meshgrid(theta_new, grid.omg_short);

EAp_new = interp2(kperp_grid, omg_grid, EAp_interp, theta_new_grid, omg_new_grid, 'linear');
EAp_zf_new = interp2(kperp_grid, omg_grid, EAp_zf_interp, theta_new_grid, omg_new_grid, 'linear');

% theta_flat = theta_grid(:);
% omg_flat = omg_grid(:);
% EAp_flat = EAp(:,i_zero+20: harm_mid_index);
% 
% F = scatteredInterpolant(theta_flat, omg_flat, EAp_flat, 'linear');
% EAp_new = F(theta_new_grid, omg_new_grid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
orders = 1e-3;

figure; 
subplot(2,2,1);
max_EAp = max(max(abs(EAp)));
imagesc(grid.omg(i_zero+20:harm_mid_index), grid.kperp, abs(EAp(:,i_zero+20:harm_mid_index))); 
axis xy; 
colorbar; colormap(parula);
colorbar; set(gca, 'ColorScale','log'); 
set(gca, 'CLim', [orders*max_EAp,max_EAp]); % Orders of magnitude
ylabel("$k_\perp [m^{-1}]$"); xlabel("$\omega [\omega_0]$");  title("$E(z=0)$")

subplot(2,2,2);
imagesc(grid.omg(i_zero+20:harm_mid_index), theta_new, EAp_new); 
axis xy; 
colorbar; colormap(parula);
colorbar; set(gca, 'ColorScale','log'); 
set(gca, 'CLim', [orders*max_EAp,max_EAp]); % Orders of magnitude
ylabel("$\theta [rad]$"); xlabel("$\omega [\omega_0]$");  title("$E(z=0)$")

subplot(2,2,3);
% Plot E(k_perp, omg, z=0) colormap on k_perp vs omg 
imagesc(grid.omg(i_zero+20:harm_mid_index), grid.kperp, abs(EAp_zf(:,i_zero+20:harm_mid_index))); 
axis xy; 
colorbar; colormap(parula);
colorbar; set(gca, 'ColorScale','log'); 
set(gca, 'CLim', [orders*max_EAp,max_EAp]); % Orders of magnitude
ylabel("$k_\perp [m^{-1}]$"); xlabel("$\omega [\omega_0]$");  title("$E(z=z_f)$")

subplot(2,2,4);
imagesc(grid.omg(i_zero+20:harm_mid_index), theta_new, abs(EAp_zf_new)); 
axis xy; 
colorbar; colormap(parula);
colorbar; set(gca, 'ColorScale','log'); 
set(gca, 'CLim', [orders*max_EAp,max_EAp]); % Orders of magnitude
ylabel("$\theta [rad]$"); xlabel("$\omega [\omega_0]$");  title("$E(z=z_f)$")



% ------------------ OLD --------------------------------
% theta_arr = linspace(0,2*pi,length(grid.omg_rad));
% [omega_2D_plot, theta_2D] = meshgrid(grid.omg_rad, theta_arr);
% 
% k_perp = (grid.omg / const.cl) .* sin(theta_arr);
% [omega_2D_temp, k_perp_2D] = meshgrid(grid.omg, grid.kperp);
% 
% % Calculate k_z = sqrt(omg^2 * n(omg) /c^2 -kperp^2) 
% % grid.kz = sqrt(grid.omg2D.^2)
% 
% 
% temp_omg = grid.omg_rad .* nindex / const.cl;
% 
% % Define meshgrid for k_perp and omg
% [grid.kperp2D,grid.omg2D] = meshgrid(grid.kperp, temp_omg);
% grid.kz = sqrt((grid.omg2D.^2) - grid.kperp2D.^2);
% 
% % Calculate theta = atan(k_perp/k_z)
% theta = atan(grid.kperp2D./grid.kz);
% --------------------------------------------------------



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