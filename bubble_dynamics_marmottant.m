%% bubble_dynamics_marmottant.m
% Bubble dynamics of a gas-core lipid-shell microbubble under ultrasound.
% Model: modified Rayleigh-Plesset + Marmottant shell law.
%
% IMPORTANT:
% - The "rupture" result means the modeled radius reached R_rupture.
%   It does NOT prove fragmentation, payload release, or biological efficacy.
% - Fluorescent surrogate (e.g., AMC) is NOT assigned a direct effect by default.
%   If measured data show a change, update shell_elasticity or shell_viscosity.
% - N2 is represented through gas polytropic behavior only. Gas diffusion,
%   dissolution, thermal effects, and shell chemistry are not modeled.
% - All defaults are illustrative starting values, not validated for your sample.

clear;
close all;
clc;

%% 1) Physical parameters (SI units)
par.R0 = 2.0e-6;            % equilibrium bubble radius [m] (2 um)
par.rho = 1000;             % surrounding liquid density [kg/m^3]
par.mu = 1.0e-3;            % liquid dynamic viscosity [Pa*s]
par.P0 = 101325;            % ambient pressure [Pa]
par.Pvap = 2330;            % water vapor pressure near room temperature [Pa]
par.kappa_gas = 1.4;        % polytropic exponent (illustrative)
par.sigma_water = 0.072;    % surface tension of clean water [N/m]

% Marmottant shell parameters (must be calibrated for the actual formulation)
par.Rbuckling = 0.95 * par.R0;  % buckling radius [m]
par.shell_elasticity = 0.50;    % shell elasticity chi [N/m]
par.shell_viscosity = 1.0e-8;   % shell dilatational viscosity [kg/s]

% For this model, the gas pressure at R0 balances ambient pressure plus
% the equilibrium shell Laplace pressure. The initial surface tension is
% obtained from the Marmottant law at R0.
sigma0 = marmottant_sigma(par.R0, par);
par.Pg0 = par.P0 - par.Pvap + 2*sigma0/par.R0; % gas pressure at R0 [Pa]

% Rupture radius in the Marmottant law: shell tension reaches sigma_water.
par.Rrupture = par.Rbuckling * sqrt(1 + par.sigma_water/par.shell_elasticity);

fprintf('R0 = %.3f um\n', par.R0*1e6);
fprintf('Modeled rupture radius = %.3f um\n', par.Rrupture*1e6);
fprintf('Note: rupture radius depends strongly on shell parameters.\n\n');

%% 2) Single run: set ultrasound conditions here
f0 = 1.0e6;             % ultrasound frequency [Hz]
Pac = 100e3;            % peak negative pressure amplitude [Pa] (100 kPa)
nCycles = 20;           % number of cycles
nPtsPerCycle = 100;     % temporal resolution for output
tEnd = nCycles/f0;
tspan = linspace(0, tEnd, nCycles*nPtsPerCycle + 1);

% State y = [R; dR/dt]
y0 = [par.R0; 0];
opts = odeset('RelTol',1e-8,'AbsTol',[1e-12 1e-6], ...
              'MaxStep',1/(f0*50));
[t,y] = ode45(@(t,y) bubble_ode(t,y,par,Pac,f0), tspan, y0, opts);

R = y(:,1);
Rdot = y(:,2);
Rmax = max(R);
Rmin = min(R);
ruptured = any(R >= par.Rrupture);

fprintf('Single run: f = %.2f MHz, Pneg = %.0f kPa\n', f0/1e6, Pac/1e3);
fprintf('Rmax/R0 = %.3f, Rmin/R0 = %.3f\n', Rmax/par.R0, Rmin/par.R0);
fprintf('Modeled shell rupture threshold reached? %s\n\n', string(ruptured));

figure('Name','Bubble radius vs time');
plot(t*1e6, R/par.R0, 'LineWidth', 1.5); hold on;
yline(par.Rrupture/par.R0, '--r', 'Modeled rupture radius');
xlabel('Time (\mus)'); ylabel('R(t) / R_0');
title(sprintf('Bubble response: %.2f MHz, %.0f kPa',f0/1e6,Pac/1e3));
grid on;

figure('Name','Bubble radius phase plot');
plot(R*1e6, Rdot, 'LineWidth', 1.2);
xlabel('Radius R (\mum)'); ylabel('Radial velocity dR/dt (m/s)');
title('Bubble radial dynamics'); grid on;

%% 3) Sweep acoustic pressure to estimate modeled shell-rupture threshold
% This is a coarse parameter sweep, not a clinical safety threshold.
pressures = (20e3:20e3:600e3);  % [Pa]
ruptureFlag = false(size(pressures));
RmaxRatio = zeros(size(pressures));

for k = 1:numel(pressures)
    Pac_k = pressures(k);
    [tk,yk] = ode45(@(t,y) bubble_ode(t,y,par,Pac_k,f0), ...
        [0, tEnd], y0, opts);
    RmaxRatio(k) = max(yk(:,1))/par.R0;
    ruptureFlag(k) = max(yk(:,1)) >= par.Rrupture;
end

idx = find(ruptureFlag,1,'first');
if isempty(idx)
    fprintf('No modeled shell rupture in tested range (up to %.0f kPa).\n', ...
        max(pressures)/1e3);
    Pthreshold = NaN;
else
    Pthreshold = pressures(idx);
    fprintf('Coarse modeled rupture threshold: about %.0f kPa at %.2f MHz.\n', ...
        Pthreshold/1e3, f0/1e6);
    fprintf('Refine pressure step and validate parameters before interpreting this value.\n');
end

figure('Name','Pressure sweep');
plot(pressures/1e3, RmaxRatio, 'o-', 'LineWidth',1.5); hold on;
yline(par.Rrupture/par.R0, '--r', 'Modeled rupture radius');
if ~isnan(Pthreshold), xline(Pthreshold/1e3, '--k', 'First crossing'); end
xlabel('Peak negative pressure (kPa)'); ylabel('Maximum R / R_0');
title(sprintf('Pressure sweep at %.2f MHz', f0/1e6)); grid on;

figure('Name','Rupture map over pressure and frequency');
freqs = [0.5e6 0.75e6 1e6 1.5e6 2e6]; % [Hz], example values
Pscan = 20e3:20e3:600e3;
ruptureMap = false(numel(freqs),numel(Pscan));
maxMap = zeros(size(ruptureMap));

for i = 1:numel(freqs)
    f = freqs(i);
    tEnd_i = nCycles/f;
    opts_i = odeset('RelTol',1e-8,'AbsTol',[1e-12 1e-6], ...
                    'MaxStep',1/(f*50));
    for j = 1:numel(Pscan)
        [~, yy] = ode45(@(tt,yy) bubble_ode(tt,yy,par,Pscan(j),f), ...
            [0 tEnd_i], y0, opts_i);
        maxMap(i,j) = max(yy(:,1))/par.R0;
        ruptureMap(i,j) = max(yy(:,1)) >= par.Rrupture;
    end
end

imagesc(Pscan/1e3, freqs/1e6, ruptureMap);
set(gca,'YDir','normal'); xlabel('Peak negative pressure (kPa)');
ylabel('Frequency (MHz)'); title('Modeled shell-rupture map (1 = threshold reached)');
colorbar; grid on;

%% ---- Local functions ----
function dydt = bubble_ode(t,y,p,Pac,f)
    R = max(y(1), 0.1*p.R0); % numerical guard only
    Rdot = y(2);

    % Acoustic pressure is positive during compression and negative during rarefaction.
    P_acoustic = Pac * sin(2*pi*f*t);

    % Gas pressure with polytropic compression/expansion.
    Pg = p.Pg0 * (p.R0/R)^(3*p.kappa_gas);

    % Marmottant shell surface tension: buckled / elastic / ruptured.
    sigma = marmottant_sigma(R,p);

    % Modified Rayleigh-Plesset equation.
    numerator = Pg - p.P0 + p.Pvap - 2*sigma/R ...
                - 4*p.mu*Rdot/R ...
                - 4*p.shell_viscosity*Rdot/(R^2) ...
                - P_acoustic;
    Rddot = numerator/(p.rho*R) - 1.5*(Rdot^2)/R;

    dydt = [Rdot; Rddot];
end

function sigma = marmottant_sigma(R,p)
    if R <= p.Rbuckling
        sigma = 0; % buckled shell
    elseif R < p.Rrupture
        sigma = p.shell_elasticity * ((R/p.Rbuckling)^2 - 1);
        sigma = max(0,min(sigma,p.sigma_water));
    else
        sigma = p.sigma_water; % shell considered ruptured
    end
end
