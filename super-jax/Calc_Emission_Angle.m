%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prelim quantities for transformations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ell = 0; % orbital angular momentum of pulse --> Bessel of order 0
grid.rmax = max(grid.r); % Get maximum R

[grid.r, grid.kperp] = DHTarrays(grid.Nr,grid.rmax,ell);
[grid.T_R2K,grid.T_K2R] = DHTtransformMatrix(grid.Nr,grid.rmax,ell);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transformations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------------- 
% Perform Hankel transform for E(r,omg,z) --> E(k_perp,omg,z)
% -------------------------------------------------------------------------
% Efinal = EF.r(:,:,end) + 1j * EF.i(:,:,end);
Efinal = Efinal_THz_t;
Ef_trans = grid.Nt * grid.dt * ifftshift(ifft(Efinal,[],2),2);
Ef_trans = transpose(grid.T_R2K*Efinal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate angle of emmission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~, i_zero] = min(abs(grid.omg));
grid.omgpos = grid.omg(i_zero:end);

[kperp_grid, omg_grid] = meshgrid(grid.kperp, grid.omgpos);
theta_grid = real(asin(const.cl * kperp_grid ./ omg_grid));
max_theta = max(theta_grid, [], 'all');

theta_g = linspace(0,max_theta,grid.Nr);
% theta_g = linspace(-0.12,0.12,grid.Nr);

% Test function
Etest = exp(-(theta_grid - 0.05).^2/(.02)^2);

if brute
    Calc_Emission_Angle_BRUTE;
else
    Calc_Emission_Angle_INTERP;
end




% Ef angle plot -----------------------------------------------------------
% figure; imagesc(grid.t, grid.r, real(Efinal)'); axis xy;
% c = redblueTecplot(); colormap(c); colorbar; set(gca, 'CLim', [-max_THz, max_THz]);
% xlabel("t [rs]"); ylabel("r [m]"); title("Efinal - Before")

% Account for dispersion?
% Keeps getting error for interp2 "Sample points must be sorted in
% ascending order."
% [kperp_grid, omg_grid] = meshgrid(grid.kperp, omg_dispersion);
% theta_grid = asin(const.cl * kperp_grid ./ (omg_grid .* nindex'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate Cherenkov Angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[max_val, max_idx] = max(real(Eg(:)));
[i_omega, i_theta] = ind2sub(size(Eg), max_idx);
theta_C = theta_g(i_theta);


% Testing
% v_source = const.cl * 1.0; % source velocity in m/s
% cos_theta_C = const.cl ./ (nindex * v_source);
% cos_theta_C = min(max(cos_theta_C, -1), 1);  % avoid domain errors
% theta_C = acos(cos_theta_C);

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