function [Signal_time] = SignalModelR2(D_axons,V_Fiber,G_axons,params)

% This function estimates the averaged signal in an unitary voxel with a
% specific amount of axons aiming to specific directions (D_axons), with 
% specific volume (V_axons) and G ratio (G_axons) dispersed with a specific
% Watson distribution coefficient (kappa) and with the main orientation 
% aligned with a specific angle (theta) with respect B0.
% A FiberLimit variable was added to define the total fibre volume in the
% voxel. Therefore, a crossing check is done: sum(V_Fiber(i),i = N_Fibers)
% = FiberLimit. Otherwise, if len(V_Fiber) = 1, then V_Fiber = 
% FiberLimit/N_Fibers. If sum(V_Fiber(i),i) > N_Fibers, then V_Fiber(i) =
% V_Fiber(i)/sum(V_Fiber(i),i) * FiberLimit/N_Fiber (keeps the same weight
% but normalised with respect the fiber limit volume).
% Most of the parameters are already set-up as for example susceptibilities
% for each compartment and total volume fractions. The model used for each
% compartment is decided below.
% Update: 11.10.2019 SNR coefficient added.
% Update: 09.12.2019 SNR removed (since it is done outside this code),
% sequence cleaning and trying to make alpha parameter free.
% Update: 05.01.2020 Myelin option removed. This code gives signal
% functions with and without myelin contribution. It avoids re-calculation
% redundancies.
% Update: 18.03.2020 Dephase_E is defined analytically from Jablinsky.

%   Detailed explanation goes here
% Input:
% *D_axons*: matrix with dimension (m,3).
% This variable contains the coordinate-vectors of all the cylinders used 
% in the simulation. These vectors are NORMALISED and in cartesian
% coordinates (x,y,z).
% *V_Fiber*: column array with dimension (m,1) or unique double value
% This variable defines the fibre volume of each cylinder (if it has size
% (m,1)) or an unique value if all the cylinders have the same volume fraction.
% *G_axons*: column array with dimension (m,1) or unique double value
% This variable defines the g-ratio value of each cylinder (if it has size
% (m,1)) or an unique value if all the cylinders have the same g-ratio.
% *params*: structure variable with the following fields listed below.
% If a field is missing then it is replaced with predefined values (see 
% next section of the code). However some fields are mandatory.
%   (MANDATORY) params.theta: mean angular orientation (in radians)
%   (MANDATORY) params.kappa: dispersion value (positive arbitrary value)
%   params.chiI: isotropic susceptibility (in ppm)
%   params.chiA: anisotropic susceptibility (in ppm)
%   params.E: exchange factor(in ppm)
%   params.B0: main magnetic field (in T)
%   params.ProtM: proton density in myelin compartment.
%   params.ProtA: proton density in intracellular compartment.
%   params.ProtE: proton density in extracellular compartment.
%   params.R2m: relaxation value in myelin compartment (in 1/ms).
%   params.R2e: relaxation value in extracellular compartment (in 1/ms).
%   params.R2a: relaxation value in intracellular compartment (in 1/ms).
%   params.FiberLimit: check that the FVF per "cylinder" or the total FVF 
%                      is not larger than this limit.
%                      
% Output:
% *Signal_time* - cell of function handles: Signal_time{1,4} @(t) 
% This is a cell of function handles that simulates the signal decay in
% function of time following the Wharton & Bowtell signal modelling per
% compartment. The following signal equations are defined:
%   Signal_time{1} @(t): Signal decay with dispersion and no myelin.
%   Signal_time{2} @(t): Signal decay with dispersion and myelin.
%   Signal_time{3} @(t): Signal decay full parallel and no myelin.
%   Signal_time{4} @(t): Signal decay full parallel and with myelin.

% References:
% 1.- Yablonskiy, D.A. and Haacke, E.M. (1994), Theory of NMR signal behavior 
% in magnetically inhomogeneous tissues: The static dephasing regime. Magn. 
% Reson. Med., 32: 749-763. doi:10.1002/mrm.1910320610
% 2.- Wharton, S. and Bowtell, R. (2012), Fiber orientation-dependent WM 
% contrast in GE MRI. PNAS, 109 (45): 18559-18564; doi: 10.1073/pnas.1211075109
% 3.- Wharton, S. and Bowtell, R. (2013), Gradient echo based fiber 
% orientation mapping using R2* and frequency difference measurements. 
% Neuroimage, 83: 1011-1023. doi:10.1016/j.neuroimage.2013.07.054

%% Library
% Here all the custom functions required by the simulation are loaded

% Struve function, used for the analytical solution of Yablonsky's Integral.
addpath('Struve01'); % This folder must be contained in the same directory that this code.

%% Internal variables
% If params variable is not defined, then the following values are defined.

if ~isfield(params,'theta') || ~isfield(params,'kappa')
    error('Mean angular orientation (params.theta) or dispersion value (params.kappa) is missing');
end

if ~isfield(params,'chiI')
    params.chiI = -0.1*10^(-6);
end

if ~isfield(params,'chiA')
    params.chiA = -0.1*10^(-6);
end

if ~isfield(params,'E')
    params.E = 0.02*10^(-6);
end

if ~isfield(params,'B0')
    params.B0 = 7.0;
end

if ~isfield(params,'ProtM')
    params.ProtM = 0.7*5000;
end

if ~isfield(params,'ProtA')
    params.ProtA = 1.0*5000;
end

if ~isfield(params,'ProtE')
    params.ProtE = 1.0*5000;
end

if ~isfield(params,'R2a') % in 1/ms
    params.R2a = 1/(0.036*1000);
end

if ~isfield(params,'R2e') % in 1/ms
    params.R2e = 1/(0.036*1000);
end

if ~isfield(params,'R2m') % in 1/ms
    params.R2m = 1/(0.008*1000);
end

if ~isfield(params,'FiberLimit')
    params.FiberLimit = 0.5;
end

% Fixed parameters
gamma = (267.52*10^(6))*1e-3; % gyromagnetic ratio in rad/(T*msec)    

%% Functions

% Auxiliar functions
% Susceptibility of the extra-cellular space (from W&B 2013)
chi_D = @(G) (params.chiI + params.chiA/4)*(1 - G.^2);
% Auxiliar function for the frequency offset in myelin compartment (from
% W&B 2013)
c_1 = @(G) 1/4 + (3.*G.^2.*log(G)).*(2.*(1-G.^2)).^(-1);
% Watson distribution: k is kappa, (mean_theta,mean_phi) are the
% spherical angles of the mean orientation and (angle_theta,angle_phi) of
% each cylinder.
Watson = @(k, mean_theta, angle_theta, mean_phi, angle_phi) (hypergeom(1/2,3/2,k)).^(-1) .* exp(k .* (sin(mean_theta).*sin(angle_theta).*cos(mean_phi - angle_phi) + cos(mean_theta).*cos(angle_theta)).^2);

% Frequency and dephase functions (from W&B 2013)
Freq_M = @(angle, G) gamma.*params.B0.*(params.chiI./2.*(2/3 - sin(angle).^2) + params.chiA./2.*(c_1(G).*sin(angle).^2 - 1/3) + params.E);
Freq_A = @(angle, G) gamma.*params.B0.*(3.*params.chiA./4.*sin(angle).^2.*log(1./G));
Freq_E = @(angle, G) gamma.*params.B0.*0.5.*abs(chi_D(G)).*sin(angle).^2;
%Int_fun = @(u,angle,G,t) (1 - besselj(0,Freq_E(angle,G).*t.*u))./u.^2;
%Int_fun, defined in Jablonskiy 1994, can be analytically defined by the
%next function (frintegral(@(u)Int_fun(u,angle,G,t),0,1)om Mathematica):
Integrated_Fun = @(angle, G, t) 0.5.*(-2 + Freq_E(angle,G).*t.*besselj(1,Freq_E(angle,G).*t).*(-2 + pi.*(Freq_E(angle,G).*t).*StruveH0(Freq_E(angle,G).*t))...
                + besselj(0,Freq_E(angle,G).*t).*(2 + 2.*(Freq_E(angle,G).*t).^2 - pi.*(Freq_E(angle,G).*t).^2.*StruveH1(Freq_E(angle,G).*t)));
Dephase_E = @(V_Fib, angle, G, t) V_Fib.* Integrated_Fun(angle,G,t);
Dephase_M = 0;

% Signal functions
S_E = @(VE, rhoE, angle, R2E, G, t) VE.*rhoE.*exp(-t.*R2E - Dephase_E(1-VE,angle,G,t));
S_A = @(VA, rhoA, angle, R2A, G, t) VA.*rhoA.*exp(-t.*R2A + 1i.*Freq_A(angle,G).*t);
S_M = @(VM, rhoM, angle, R2M, G, t) VM.*rhoM.*exp(-t.*R2M + 1i.*Freq_M(angle,G).*t - Dephase_M);

%% Input variables
 
[phi_axons, theta_axons, ~] = cart2sph(D_axons(:,1), D_axons(:,2), D_axons(:,3));
theta_axons = pi/2 - theta_axons; % Matlab defines pi/2 in the +z axis.

if numel(V_Fiber) == 1
    V_axons = V_Fiber*ones(numel(D_axons(:,1)),1);
else
    V_axons = V_Fiber;
    if any(V_Fiber > FiberLimit)
        V_axons(V_Fiber > FiberLimit) = FiberLimit;
    end
end

% The Watson distribution is evaluated by the cylinder's spherical angles
% it is assumed that the azimutal angle of the mean orientation is 0.
Weight_Watson = @(k, mean_theta) Watson(k, mean_theta, theta_axons, 0, phi_axons);
Watson_Norm_Factor = @(k, mean_theta) sum(Weight_Watson(k, mean_theta));

if numel(G_axons) == 1
    G_axons = G_axons * ones(numel(D_axons(:,1)),1);
else
    if max(G_axons) > 1.0
        G_axons = G_axons/max(G_axons);
    elseif max(G_axons) < 0.0
        G_axons = (G_axons + min(G_axons))/(max(G_axons) - min(G_axons));
    end
end

V_axons_M = V_axons.*(1 - G_axons.^2); % Myelin volume per axon
V_axons_A = V_axons.*(G_axons.^2); % Axonal volume per axon
V_extra = 1 - V_axons; % Extracellular volume is the sum of the rest

%% Signal calculation

% Estimate the contribution of each compartment at each theta_r and phi_r
% and kappa:

S_E_kappa = @(angle, kdis, t) sum(S_E(V_extra, params.ProtE, theta_axons, params.R2e, G_axons, t).*Weight_Watson(kdis,angle));
S_A_kappa = @(angle, kdis, t) sum(S_A(V_axons_A, params.ProtA, theta_axons, params.R2a, G_axons, t).*Weight_Watson(kdis,angle));
S_M_kappa = @(angle, kdis, t) sum(S_M(V_axons_M, params.ProtM, theta_axons, params.R2m, G_axons, t).*Weight_Watson(kdis,angle));

Signal_time{1} = @(t) (S_E_kappa(params.theta, params.kappa, t) + ...
                       S_A_kappa(params.theta, params.kappa, t))./Watson_Norm_Factor(params.kappa, params.theta);

Signal_time{2} = @(t) (S_E_kappa(params.theta, params.kappa, t) + ...
                       S_A_kappa(params.theta, params.kappa, t) + ...
                       S_M_kappa(params.theta, params.kappa, t))./Watson_Norm_Factor(params.kappa, params.theta);

% Parallel case, where all the cylinders are parallel to the main
% orientation. However, with "all" the cylinders I refer to equal number of
% cylinders used in the dispersion case.
theta_parallel = params.theta * ones(numel(D_axons(:,1)),1);

S_E_no_kappa = @(t) sum(S_E(V_extra, params.ProtE, theta_parallel, params.R2e, G_axons, t));
S_A_no_kappa = @(t) sum(S_A(V_axons_A, params.ProtA, theta_parallel, params.R2a, G_axons, t));
S_M_no_kappa = @(t) sum(S_M(V_axons_M, params.ProtM, theta_parallel, params.R2m, G_axons, t));

% The number of cylinders must be the same than the used in the dispersion
% case!
Signal_time{3} = @(t) (S_E_no_kappa(t) + ...
                       S_A_no_kappa(t))./numel(theta_parallel);
                          
Signal_time{4} = @(t) (S_E_no_kappa(t) + ...
                       S_A_no_kappa(t) + ...
                       S_M_no_kappa(t))./numel(theta_parallel);
end