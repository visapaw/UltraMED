%  MICROBU BBLE ULTRASOUND SIMULATION
%  Lipid shell + AMC fluorescent surrogate + N2 gas core
%  Assumption using
%   1. microbubbleมีรูปร่างเป็นทรงกลมเสมอ
%   2. ก๊าซภายในฟองไม่มีการรั่วไหลจนถึงRupture
%   3. ของเหลวมีคุณสมบัติคงที่ เช่น ความหนาแน่นและความหนืด
%   4. เปลือกไขมันมีสมบัติเชิงกลคงที่ และใช้แบบจำลอง Marmottant-type อธิบายแรงตึงผิวของเปลือก
%   5. คลื่นอัลตราซาวด์เป็นคลื่นไซน์ ที่มีความถี่และแอมพลิจูดความดันคงที่
%   6. ไม่พิจารณาพฤติกรรมหลังmicrobubbleแตก
%   7. ไม่พิจารณาปฏิสัมพันธ์ระหว่างฟอง และผลจากการเปลี่ยนแปลงอุณหภูมิ

clear; clc; close all;

%% ============================================================
%  SECTION 1: USER-DEFINED / OPTIMIZATION PARAMETERS

DPPC.mg = 30.05;
DPPA.mg = 5.00;
DPPE.mg = 10.00;
DSPE.mg = 50.05;
lipid.DPPC = DPPC.mg/(DPPC.mg+DPPA.mg+DPPE.mg+DSPE.mg);
lipid.DPPA = DPPA.mg/(DPPC.mg+DPPA.mg+DPPE.mg+DSPE.mg);
lipid.DPPE = DPPE.mg/(DPPC.mg+DPPA.mg+DPPE.mg+DSPE.mg);
lipid.DSPE = DSPE.mg/(DPPC.mg+DPPA.mg+DPPE.mg+DSPE.mg);

lipidRatio = [lipid.DPPC, lipid.DPPA, ...
              lipid.DPPE, lipid.DSPE];
assert(all(lipidRatio >= 0), ...
    'Lipid ratios must be non-negative.');
assert(abs(sum(lipidRatio) - 1) < 1e-8, ...
    'Lipid ratios must sum to 1.');
disp(lipidRatio);

% ---- Ultrasound parameters ----
US.f       = 0.25e6;       % Hz
US.Pa      = 150000;       % Pa; peak acoustic pressure amplitude
US.nCycles = 30;          % Number of ultrasound cycles
US.phase   = 0;           % rad

% ---- Shell mechanical parameters ----
% the specified DPPC/DPPA/DPPE/DSPE + AMC formulation.
p.chi = 0.53;           % N/m; shell elasticity
p.kappaS = 1e-8;       % N*s/m; shell dilatational viscosity
p.sigma0 = 0;

% Rupture threshold:
% Set to NaN to calculate a Marmottant-style estimate.
shell.Rrupture = NaN;      % m
% Optional measured buckling radius.
shell.Rbuckling = NaN;     % m

% ---- Gas model ----
gas.gamma = 1.4;           % Approximate polytropic exponent for N2

% AMC is treated as a fluorescent surrogate/tracer here.
% It does not automatically change gas pressure or shell mechanics.
AMC.enabled = true;        % Metadata only; no photophysics modeled

%% ============================================================

%  SECTION 2: FIXED PHYSICAL PARAMETERS
R0 = 1e-6;               % m
p.P0 = 1.046e5;             % Pa; ambient pressure
p.rho = 1000;               % kg/m^3; liquid density
p.mu = 0.004;             % Pa*s; liquid dynamic viscosity
p.c = 1500;                % m/s; speed of sound in liquid
p.sigmaWater = 0.072;      % N/m
p.gasGamma = gas.gamma;
p.R0 = R0;

% Buckling radius
if isnan(shell.Rbuckling)
    p.Rbuckling = R0 / sqrt(1 + p.sigma0 / p.chi);
else
    p.Rbuckling = shell.Rbuckling;
end

% Rupture radius
if isnan(shell.Rrupture)
    p.Rrupture = p.Rbuckling * ...
        sqrt(1 + p.sigmaWater / p.chi);
else
    p.Rrupture = shell.Rrupture;
end

% Check model configuration
assert(R0 > 0, 'R0 must be positive.');
assert(US.f > 0, 'Frequency must be positive.');
assert(US.Pa >= 0, 'Pressure amplitude cannot be negative.');
assert(p.chi > 0, 'Shell elasticity must be positive.');
assert(p.kappaS >= 0, 'Shell viscosity cannot be negative.');
assert(p.Rrupture > p.Rbuckling, ...
    'Rupture radius must exceed buckling radius.');
%% ============================================================


%% ============================================================
% SECTION 3: SIMULATE BOTH RP AND KELLER-MIKSIS

T = 1 / US.f;
tEnd = US.nCycles * T;

y0 = [R0; 0];
tOutput = linspace(0, tEnd, 4000);

options = odeset('RelTol', 1e-6, ...
                 'AbsTol', [1e-11, 1e-4], ...
                 'MaxStep', T/200);

% ---- RP simulation ----
[tRP, yRP] = ode45( ...
    @(t,y) bubbleODE_RP(t, y, p, US), ...
    tOutput, y0, options);

RRP = yRP(:,1);
URP = yRP(:,2);

% ---- Keller-Miksis simulation ----
[tKM, yKM] = ode45( ...
    @(t,y) bubbleODE_KM(t, y, p, US), ...
    tOutput, y0, options);

RKM = yKM(:,1);
UKM = yKM(:,2);

% ---- Gilmore simulation ----
[tGil, yGil] = ode15s( ...
    @(t,y) bubbleODE_Gilmore(t, y, p, US), ...
    tOutput, y0, options);

RGil = yGil(:,1);
UGil = yGil(:,2);

% ---- Calculate shell surface tension ----
sigmaRP  = zeros(size(RRP));
sigmaKM  = zeros(size(RKM));
sigmaGil = zeros(size(RGil));

for k = 1:length(RRP)
    sigmaRP(k) = shellModel(RRP(k), p, false);
end

for k = 1:length(RKM)
    sigmaKM(k) = shellModel(RKM(k), p, false);
end

for k = 1:length(RGil)
    sigmaGil(k) = shellModel(RGil(k), p, false);
end

fprintf('\n========== RP vs KM vs GILMORE ==========\n');
fprintf('RP maximum radius      = %.4f um\n', max(RRP)*1e6);
fprintf('KM maximum radius      = %.4f um\n', max(RKM)*1e6);
fprintf('Gilmore maximum radius = %.4f um\n', max(RGil)*1e6);

fprintf('RP minimum radius      = %.4f um\n', min(RRP)*1e6);
fprintf('KM minimum radius      = %.4f um\n', min(RKM)*1e6);
fprintf('Gilmore minimum radius = %.4f um\n', min(RGil)*1e6);

fprintf('RP reaches threshold      = %d\n', any(RRP >= p.Rrupture));
fprintf('KM reaches threshold      = %d\n', any(RKM >= p.Rrupture));
fprintf('Gilmore reaches threshold = %d\n', any(RGil >= p.Rrupture));

% Colors
colorRP  = [0 0.75 0.85];
colorKM  = [0.95 0.45 0.10];
colorGil = [0.55 0.30 0.75];

%% ============================================================
% PLOT 1: RADIUS VS TIME

figure('Color','w','Name','RP vs KM vs Gilmore - Radius');
hold on; box on; grid on;

plot(tRP*1e6, RRP*1e6, ...
    'Color', colorRP, 'LineWidth', 1.5);

plot(tKM*1e6, RKM*1e6, ...
    'Color', colorKM, 'LineWidth', 1.5);

plot(tGil*1e6, RGil*1e6, ...
    'Color', colorGil, 'LineWidth', 1.5);

yline(R0*1e6, '--k', 'R_0');
yline(p.Rrupture*1e6, '--r', 'Rupture threshold');

xlabel('Time (\mus)');
ylabel('Bubble radius R (\mum)');
title('Bubble Oscillation: RP vs KM vs Gilmore');

legend('RP + Marmottant', ...
       'KM + Marmottant', ...
       'Gilmore + Marmottant', ...
       'Initial radius', ...
       'Rupture threshold', ...
       'Location', 'best');

%% ============================================================
% PLOT 2: SHELL SURFACE TENSION VS RADIUS


figure('Color','w','Name','RP vs KM vs Gilmore - Shell');
hold on; box on; grid on;

plot(RRP*1e6, sigmaRP*1e3, ...
    'Color', colorRP, 'LineWidth', 1.4);

plot(RKM*1e6, sigmaKM*1e3, ...
    'Color', colorKM, 'LineWidth', 1.4);

plot(RGil*1e6, sigmaGil*1e3, ...
    'Color', colorGil, 'LineWidth', 1.4);

xline(p.Rbuckling*1e6, '--b', 'Buckling radius');
xline(p.Rrupture*1e6, '--r', 'Rupture threshold');

xlabel('Bubble radius R (\mum)');
ylabel('Effective surface tension (mN/m)');
title('Marmottant Shell Response: RP vs KM vs Gilmore');

legend('RP + Marmottant', ...
       'KM + Marmottant', ...
       'Gilmore + Marmottant', ...
       'Buckling radius', ...
       'Rupture threshold', ...
       'Location', 'best');

%% ============================================================
% PLOT 3: PRESSURE AMPLITUDE SWEEP

PaList = linspace(10e3, 500e3, 300);

maxR_RP = zeros(size(PaList));
maxR_KM = zeros(size(PaList));
maxR_Gil = zeros(size(PaList));

ruptureRP = false(size(PaList));
ruptureKM = false(size(PaList));
ruptureGil = false(size(PaList));

for j = 1:length(PaList)

    USj = US;
    USj.Pa = PaList(j);

    % RP
    [~, yjRP] = ode45( ...
        @(t,y) bubbleODE_RP(t, y, p, USj), ...
        tOutput, y0, options);

    maxR_RP(j) = max(yjRP(:,1));
    ruptureRP(j) = any(yjRP(:,1) >= p.Rrupture);

    % Keller-Miksis
    [~, yjKM] = ode45( ...
        @(t,y) bubbleODE_KM(t, y, p, USj), ...
        tOutput, y0, options);

    maxR_KM(j) = max(yjKM(:,1));
    ruptureKM(j) = any(yjKM(:,1) >= p.Rrupture);

    % Gilmore
    [~, yjGil] = ode15s( ...
        @(t,y) bubbleODE_Gilmore(t, y, p, USj), ...
        tOutput, y0, options);

    maxR_Gil(j) = max(yjGil(:,1));
    ruptureGil(j) = any(yjGil(:,1) >= p.Rrupture);

end

figure('Color','w','Name','RP vs KM vs Gilmore - Pressure Sweep');
hold on; box on; grid on;

plot(PaList/1e3, maxR_RP*1e6, ...
    '-', 'Color', colorRP, 'LineWidth', 1.6);

plot(PaList/1e3, maxR_KM*1e6, ...
    '-', 'Color', colorKM, 'LineWidth', 1.6);

plot(PaList/1e3, maxR_Gil*1e6, ...
    '-', 'Color', colorGil, 'LineWidth', 1.6);

yline(p.Rrupture*1e6, '--r', 'Rupture threshold');

xlabel('Ultrasound pressure amplitude (kPa)');
ylabel('Maximum bubble radius (\mum)');
title('Pressure Sweep: RP vs KM vs Gilmore');

legend('RP + Marmottant', ...
       'KM + Marmottant', ...
       'Gilmore + Marmottant', ...
       'Rupture threshold', ...
       'Location', 'best');

%% ============================================================
%  SECTION 5: LOCAL FUNCTIONS


function dydt = bubbleODE_RP(t, y, p, US)

    R = max(y(1), 0.05*p.R0);
    U = y(2);

    Pac = US.Pa*sin(2*pi*US.f*t + US.phase);

    [sigmaR, ~] = shellModel(R, p, false);

    Pg = (p.P0 + 2*p.sigma0/p.R0) ...
         * (p.R0/R)^(3*p.gasGamma);

    acceleration = ( ...
        Pg - p.P0 - Pac ...
        - 2*sigmaR/R ...
        - 4*p.mu*U/R ...
        - 4*p.kappaS*U/R^2 ...
        )/(p.rho*R) ...
        - 1.5*U^2/R;

    dydt = [U; acceleration];

end


function dydt = bubbleODE_KM(t, y, p, US)

    R = max(y(1), 0.05*p.R0);
    U = y(2);

    omega = 2*pi*US.f;

    Pac = US.Pa*sin(omega*t + US.phase);
    dPacdt = US.Pa*omega*cos(omega*t + US.phase);

    [sigmaR, ~] = shellModel(R, p, false);

    % d(sigma)/dR for the piecewise Marmottant model
    if R <= p.Rbuckling || R >= p.Rrupture
        dsigma_dR = 0;
    else
        dsigma_dR = 2*p.chi*R/p.Rbuckling^2;
    end

    % Polytropic gas pressure
    % No extra (1 - 3*gamma*U/c) factor here.
    Pg = (p.P0 + 2*p.sigma0/p.R0) ...
         * (p.R0/R)^(3*p.gasGamma);

    % Net pressure at the bubble wall
    G = Pg ...
        - 2*sigmaR/R ...
        - 4*p.mu*U/R ...
        - 4*p.kappaS*U/R^2 ...
        - p.P0 - Pac;

    % Partial derivative of G with respect to R
    dGdR = -3*p.gasGamma*Pg/R ...
           - 2*dsigma_dR/R ...
           + 2*sigmaR/R^2 ...
           + 4*p.mu*U/R^2 ...
           + 8*p.kappaS*U/R^3;

    % Partial derivative of G with respect to U
    dGdU = -4*p.mu/R - 4*p.kappaS/R^2;

    % Explicit time derivative of G
    dGdt = -dPacdt;

    % Keller-Miksis equation, solved for acceleration
    numerator = (1 + U/p.c)*G ...
        + (R/p.c)*(dGdR*U + dGdt) ...
        - 1.5*p.rho*U^2*(1 - U/(3*p.c));

    denominator = p.rho*R*(1 - U/p.c) ...
        - (R/p.c)*dGdU;

    acceleration = numerator/denominator;

    dydt = [U; acceleration];

end


function dydt = bubbleODE_Gilmore(t, y, p, US)

    R = max(y(1), 0.05*p.R0);
    U = y(2);

    omega = 2*pi*US.f;
    Pac = US.Pa*sin(omega*t + US.phase);

    % Tait EOS parameters for water
    B = 304.6e6;       % Pa
    nTait = 7.0;

    % Surface tension from Marmottant model
    sigmaR = shellModel(R, p, false);

    % Polytropic gas pressure
    Pg = (p.P0 + 2*p.sigma0/p.R0) ...
         * (p.R0/R)^(3*p.gasGamma);

    % Liquid pressure at bubble wall
    pL = Pg ...
         - 2*sigmaR/R ...
         - 4*p.mu*U/R ...
         - 4*p.kappaS*U/R^2;

    % Far-field absolute pressure
    pInf = p.P0 + Pac;

    % Enthalpy difference and local sound speed
    [H, C] = gilmoreHC(pL, pInf, p.rho, p.P0, B, nTait);

    % Numerical partial derivatives of H
    dR = max(1e-4*R, 1e-11);
    dU = 1e-4;

    Rminus = max(R-dR, 0.05*p.R0);
    Rplus = R+dR;

    HR = (gilmoreHonly(Rplus,U,t,p,US,B,nTait) ...
        - gilmoreHonly(Rminus,U,t,p,US,B,nTait)) ...
        /(Rplus-Rminus);

    HU = (gilmoreHonly(R,U+dU,t,p,US,B,nTait) ...
        - gilmoreHonly(R,U-dU,t,p,US,B,nTait)) ...
        /(2*dU);

    dt = max(1e-5/US.f, 1e-10);

    Ht = (gilmoreHonly(R,U,t+dt,p,US,B,nTait) ...
        - gilmoreHonly(R,U,t-dt,p,US,B,nTait)) ...
        /(2*dt);

    % Gilmore equation
    A = R*(1-U/C);
    Q = (R/C)*(1-U/C);

    numerator = H*(1+U/C) ...
        + Q*(HR*U + Ht) ...
        - 1.5*U^2*(1-U/(3*C));

    denominator = A - Q*HU;

    if ~isfinite(denominator) || abs(denominator) < 1e-12
        error('Gilmore equation became numerically singular.');
    end

    acceleration = numerator/denominator;

    dydt = [U; acceleration];

end


function H = gilmoreHonly(R,U,t,p,US,B,nTait)

    R = max(R, 0.05*p.R0);

    Pac = US.Pa*sin(2*pi*US.f*t + US.phase);

    sigmaR = shellModel(R,p,false);

    Pg = (p.P0 + 2*p.sigma0/p.R0) ...
         * (p.R0/R)^(3*p.gasGamma);

    pL = Pg ...
         - 2*sigmaR/R ...
         - 4*p.mu*U/R ...
         - 4*p.kappaS*U/R^2;

    pInf = p.P0 + Pac;

    [H,~] = gilmoreHC(pL,pInf,p.rho,p.P0,B,nTait);

end

function [H,C] = gilmoreHC(pL,pInf,rho0,P0,B,nTait)

    % Tait equation of state
    if pL+B <= 0 || pInf+B <= 0
        error('Pressure is outside the valid Tait EOS range.');
    end

    rhoInf = rho0*((pInf+B)/(P0+B))^(1/nTait);
    cInf = sqrt(nTait*(pInf+B)/rhoInf);

    ratio = (pL+B)/(pInf+B);

    H = nTait/(nTait-1) ...
        * (pInf+B)/rhoInf ...
        * (ratio^((nTait-1)/nTait)-1);

    C = sqrt(max(cInf^2+(nTait-1)*H,1));

end

function [sigmaR, state] = shellModel(R, p, isRuptured)

    % Marmottant-type piecewise shell law.
    % The rupture state is supplied by the caller if a persistent
    % rupture simulation is implemented.

    if isRuptured || R >= p.Rrupture
        sigmaR = p.sigmaWater;
        state = "Ruptured";

    elseif R <= p.Rbuckling
        sigmaR = 0;
        state = "Buckled";

    else
        sigmaR = p.chi * ...
            (R^2/p.Rbuckling^2 - 1);
        sigmaR = max(0, min(sigmaR, p.sigmaWater));
        state = "Elastic";
    end
end

%% ============================================================