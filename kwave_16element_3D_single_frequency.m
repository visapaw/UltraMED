% AXES USED IN THIS SCRIPT
%   x = propagation direction, measured from the transducer face [mm]
%   y = lateral direction across the 16-element array [mm]
%   z = elevation / element-height direction [mm]
%
% PLOTS
%   Figure 1: peak absolute pressure on the central elevation plane (x-y)
%   Figure 2: peak negative pressure (PNP) on the same x-y plane
%   Figure 3: pressure profiles on the array centre line (y ~= 0, z ~= 0)
%
% IMPORTANT MODEL LIMITATIONS
%   1) This is a 3D acoustic model of a prescribed transducer source.
%      It does not solve piezoelectric electromechanical coupling.
%   2) sourcePressureNominal_Pa is an illustrative source-strength scale,
%      converted to approximate particle velocity using p/(rho*c). It is
%      NOT a calibrated prediction of real transducer face pressure.
%   3) The medium is homogeneous, lossless, and linear in this starter model.
%   4) The pressure maps are 2D slices from a 3D simulation, at z ~= 0.
%   5) Grid sizes are a starting point. Check memory and convergence before
%      using the results for design decisions.

clear;
close all;
clc;

%% TRANSDUCER PARAMETERS

c0   = 1500;             % Sound speed [m/s]
rho0 = 1000;             % Density [kg/m^3]
f0 = 0.25e6;              % Centre frequency [Hz] (1 MHz)
nCycles = 30;             % Number of cycles in the tone burst [-]

% Array geometry (physical units)
nElements = 16;           % Number of elements [-]
elementWidth_m  = 0.40e-3; % Element width across array, along y [m]
elementHeight_m = 3.00e-3; % Element height/elevation, along z [m]
kerf_m = 0.4e-3;         % Gap between adjacent elements [m]
pitch_m = elementWidth_m + kerf_m; % Element-centre spacing [m]

% Beamforming / target location
focusDepth_m = 12e-3;     % Electronic focus distance from face along +x [m]
steeringAngle_deg = 0;    % 0 = beam along +x; nonzero steers laterally [deg]
sourcePressureNominal_Pa = 200e3; % Illustrative nominal pressure scale [Pa]

%% NUMERICAL GRID AND TIME SETTINGS

% Grid spacing: 0.10 mm gives 15 points per wavelength at 1 MHz in water.
dx = 0.10e-3;             % Grid spacing along propagation x [m]
dy = 0.10e-3;             % Grid spacing across array y [m]
dz = 0.10e-3;             % Grid spacing along elevation z [m]

% Computational grid dimensions (number of grid points)
% Increase these only after the starter case runs successfully.
Nx = 180;                 % x-direction: propagation
Ny = 180;                 % y-direction: lateral array dimension
Nz = 72;                  % z-direction: elevation dimension

PMLSize = 10;             % PML thickness in grid points on each side
CFL = 0.15;               % CFL number used to choose time step [-]
tEnd = 22e-6;             % Total simulation time [s]

% Basic derived quantities
lambda_m = c0 / f0;        % Wavelength in medium [m]
dt = CFL * min([dx, dy, dz]) / c0;
Nt = round(tEnd / dt) + 1;

kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
kgrid.setTime(Nt, dt);

medium.sound_speed = c0;
medium.density = rho0;
% No attenuation or nonlinear parameter is defined in this starter model.

%% CONVERT ELEMENT GEOMETRY FROM METRES TO GRID POINTS

elementWidthPts  = max(1, round(elementWidth_m / dy));
elementHeightPts = max(1, round(elementHeight_m / dz));
kerfPts           = max(0, round(kerf_m / dy));
aperturePts = nElements * elementWidthPts + (nElements - 1) * kerfPts;

% Actual represented dimensions after grid rounding (for reporting)
elementWidth_grid_m  = elementWidthPts * dy;
elementHeight_grid_m = elementHeightPts * dz;
kerf_grid_m          = kerfPts * dy;
pitch_grid_m         = (elementWidthPts + kerfPts) * dy;
aperture_grid_m      = aperturePts * dy;

assert(Nx > 2*PMLSize && Ny > 2*PMLSize && Nz > 2*PMLSize, ...
    'Grid is too small for the selected PMLSize.');
assert(aperturePts <= Ny - 2*PMLSize, ...
    'Array aperture does not fit between the lateral PML regions.');
assert(elementHeightPts <= Nz - 2*PMLSize, ...
    'Element height does not fit between the elevation PML regions.');

%% DEFINE THE 3D TRANSDUCER SOURCE

xSource = PMLSize + 3;
yStart = round(Ny/2 - aperturePts/2);
zStart = round(Nz/2 - elementHeightPts/2);

assert(yStart >= PMLSize + 1 && ...
       yStart + aperturePts - 1 <= Ny - PMLSize, ...
       'Array touches or overlaps the lateral PML; adjust grid or geometry.');
assert(zStart >= PMLSize + 1 && ...
       zStart + elementHeightPts - 1 <= Nz - PMLSize, ...
       'Array touches or overlaps the elevation PML; adjust grid or geometry.');

% Create a finite-duration tone burst. kWaveTransducer uses particle velocity
% as its input signal, so scale the dimensionless tone burst by p/(rho*c).
inputSignal = toneBurst(1/dt, f0, nCycles);
inputSignal = inputSignal(:).';
inputSignal = (sourcePressureNominal_Pa / (rho0 * c0)) .* inputSignal;

transducerInput.number_elements = nElements;
transducerInput.element_width = elementWidthPts;
transducerInput.element_length = elementHeightPts;
transducerInput.element_spacing = kerfPts;
transducerInput.position = [xSource, yStart, zStart];
transducerInput.radius = inf;                 % Flat linear array
transducerInput.sound_speed = c0;
transducerInput.focus_distance = focusDepth_m;
transducerInput.steering_angle = steeringAngle_deg;
transducerInput.transmit_apodization = 'Rectangular';
transducerInput.input_signal = inputSignal;

% Construct the 3D k-Wave transducer object.
transducer = kWaveTransducer(kgrid, transducerInput);

% Centre of the aperture in grid coordinates (used for plotting profiles).
arrayCenterY = yStart + (aperturePts - 1)/2;
centralYIdx = round(arrayCenterY);
centralZIdx = round(Nz/2);

%% SENSOR: RECORD A CENTRAL ELEVATION SLICE ONLY

% The calculation is 3D, but recording only this plane saves memory.
% xRange excludes the source face and the far-side PML.
xRange = (xSource + 1):(Nx - PMLSize);
yRange = (PMLSize + 1):(Ny - PMLSize);
sensor.mask = false(Nx, Ny, Nz);
sensor.mask(xRange, yRange, centralZIdx) = true;
sensor.record = {'p_max', 'p_min', 'p'};

input_args = { ...
    'PMLSize', PMLSize, ...
    'PMLInside', true, ...
    'PlotSim', false, ...
    'DataCast', 'single'};

%% RUN THE 3D ACOUSTIC SIMULATION

fprintf('Running 3D k-Wave simulation...\n');
fprintf('Frequency: %.2f MHz | Grid: %d x %d x %d | Nt: %d\n', ...
    f0/1e6, Nx, Ny, Nz, Nt);

sensor_data = kspaceFirstOrder3D( ...
    kgrid, medium, transducer, sensor, input_args{:});

%% RESHAPE SENSOR OUTPUT AND CALCULATE PRESSURE METRICS

% Sensor points are ordered according to MATLAB linear indexing of sensor.mask.
% Because the sensor mask covers xRange-by-y at one z index, reshape to [x,y].
pMaxMap = reshape(double(sensor_data.p_max), numel(xRange), numel(yRange));
pMinMap = reshape(double(sensor_data.p_min), numel(xRange), numel(yRange));

% Peak absolute acoustic pressure and peak negative pressure (PNP)
pAbsMap = max(pMaxMap, abs(pMinMap));
pnpMap = max(0, -pMinMap);

% Axes relative to the transducer face and array centre
xAxis_mm = (xRange - xSource) * dx * 1e3;          % Propagation distance [mm]
yAxis_mm = (yRange - arrayCenterY) * dy * 1e3;     % Lateral offset [mm]
centralYSensorIdx = find(yRange == centralYIdx, 1);

%%  FIGURE 1: PEAK ABSOLUTE PRESSURE MAP (CENTRAL ELEVATION SLICE)

figure('Name', '3D k-Wave: Peak absolute pressure, central z slice');
imagesc(xAxis_mm, yAxis_mm, (pAbsMap.'/1e3));
axis image;
set(gca, 'YDir', 'normal');
xlabel('Propagation distance from transducer face, x (mm)');
ylabel('Lateral position relative to array centre, y (mm)');
title(sprintf('Peak absolute pressure, central z slice, f = %.2f MHz', f0/1e6));
cb = colorbar;
ylabel(cb, 'Peak |p| (kPa)');
hold on;
plot(focusDepth_m*1e3, 0, 'wo', 'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', 'k', 'MarkerSize', 7);
hold off;

%%  FIGURE 2: PEAK NEGATIVE PRESSURE MAP (IMPORTANT FOR MICROBUBBLES)

figure('Name', '3D k-Wave: Peak negative pressure, central z slice');
imagesc(xAxis_mm, yAxis_mm, (pnpMap.'/1e3));
axis image;
set(gca, 'YDir', 'normal');
xlabel('Propagation distance from transducer face, x (mm)');
ylabel('Lateral position relative to array centre, y (mm)');
title(sprintf('Peak negative pressure (PNP), central z slice, f = %.2f MHz', f0/1e6));
cb = colorbar;
ylabel(cb, 'PNP = -min(p), when min(p) < 0 (kPa)');
hold on;
plot(focusDepth_m*1e3, 0, 'wo', 'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', 'k', 'MarkerSize', 7);
hold off;

%%  FIGURE 3: AXIAL PRESSURE PROFILE THROUGH THE ARRAY CENTRE

% Profile at y ~= 0 and z ~= 0. PNP is especially relevant to the MB model.
profileAbs_kPa = pAbsMap(:, centralYSensorIdx) / 1e3;
profilePNP_kPa = pnpMap(:, centralYSensorIdx) / 1e3;

figure('Name', 'Axial pressure profile at array centre');
plot(xAxis_mm, profileAbs_kPa, 'k-', 'LineWidth', 1.5);
hold on;
plot(xAxis_mm, profilePNP_kPa, 'r-', 'LineWidth', 1.5);
xline(focusDepth_m*1e3, '--b', 'Nominal focus', 'LineWidth', 1.2);
hold off;
xlabel('Propagation distance from transducer face, x (mm)');
ylabel('Pressure (kPa)');
title(sprintf('Axial profile at array centre, f = %.2f MHz', f0/1e6));
legend('Peak absolute pressure, max(|p|)', 'Peak negative pressure (PNP)', ...
    'Location', 'best');
grid on;

%% Instantaneous pressure field: interference pattern

tSnapshot = 15e-6;  % Snapshot time [s]
[~, idx] = min(abs(kgrid.t_array - tSnapshot));

% Extract pressure at the selected time
pSnapshot = reshape( ...
    double(sensor_data.p(:, idx)), ...
    numel(xRange), numel(yRange));

% Convert Pa to kPa
pSnapshot_kPa = pSnapshot / 1e3;

% Symmetric color range for positive and negative pressure
pLim = max(abs(pSnapshot_kPa(:)));

figure('Name', 'Instantaneous pressure field');
imagesc(xAxis_mm, yAxis_mm, pSnapshot_kPa.');
axis image;
set(gca, 'YDir', 'normal');

xlabel('Propagation distance from transducer face, x (mm)');
ylabel('Lateral position relative to array centre, y (mm)');
title(sprintf('Instantaneous acoustic pressure, f = %.2f MHz, t = %.1f us', ...
    f0/1e6, tSnapshot*1e6));

colormap(turbo);
clim([-pLim, pLim]);

cb = colorbar;
ylabel(cb, 'Instantaneous pressure p (kPa)');

%%  REPORT KEY VALUES NEAR THE NOMINAL FOCUS

[~, focusXIdx] = min(abs(xAxis_mm - focusDepth_m*1e3));
[~, centreYIdx] = min(abs(yAxis_mm));
PabsFocus_kPa = pAbsMap(focusXIdx, centreYIdx) / 1e3;
PNPFocus_kPa = pnpMap(focusXIdx, centreYIdx) / 1e3;

fprintf('\n--- Simulation summary ---\n');
fprintf('Frequency                         : %.3f MHz\n', f0/1e6);
fprintf('Number of elements                : %d\n', nElements);
fprintf('Requested element width           : %.3f mm\n', elementWidth_m*1e3);
fprintf('Grid-represented element width    : %.3f mm\n', elementWidth_grid_m*1e3);
fprintf('Requested element height          : %.3f mm\n', elementHeight_m*1e3);
fprintf('Grid-represented element height   : %.3f mm\n', elementHeight_grid_m*1e3);
fprintf('Requested kerf                    : %.3f mm\n', kerf_m*1e3);
fprintf('Grid-represented pitch            : %.3f mm\n', pitch_grid_m*1e3);
fprintf('Grid-represented array aperture   : %.3f mm\n', aperture_grid_m*1e3);
fprintf('Nominal focus depth               : %.2f mm\n', focusDepth_m*1e3);
fprintf('Peak |p| near nominal focus       : %.2f kPa\n', PabsFocus_kPa);
fprintf('PNP near nominal focus            : %.2f kPa\n', PNPFocus_kPa);
fprintf(['NOTE: These pressure values are for the illustrative source scaling; ', ...
         'calibrate the source before making physical design conclusions.\n']);
