%% kwave_16element_linear_array.m
% Acoustic-field simulation for a 16-element linear array using k-Wave.
%
% REQUIREMENTS:
% - MATLAB
% - k-Wave toolbox installed and added to MATLAB path
%
% IMPORTANT MODEL LIMITATION:
% This models an imposed acoustic source at each element. It does NOT solve
% piezoelectric electromechanical coupling or the vibration of soft piezoceramic.
% Use measured/calibrated element surface velocity or pressure as the source.
% For actual piezoceramic design, use piezoelectric FEM (e.g., COMSOL/PZFlex)
% and then export the surface motion/pressure into an acoustic model.
%
% Grid and element dimensions below are illustrative. Adapt to the real probe.

clear;
close all;
clc;

%% 1) Medium and simulation grid
c0 = 1500;           % sound speed in water/tissue-like medium [m/s]
rho0 = 1000;         % density [kg/m^3]
f0 = 1.0e6;          % center frequency [Hz]
sourcePressure = 100e3; % nominal source pressure amplitude [Pa], illustrative
nCycles = 8;         % source burst cycles

dx = 20e-6;          % grid spacing x [m]
dy = 20e-6;          % grid spacing y [m]
Nx = 512;            % lateral grid points
Ny = 768;            % depth grid points

kgrid = kWaveGrid(Nx, dx, Ny, dy);
medium.sound_speed = c0;
medium.density = rho0;

% CFL-based time step; k-Wave uses the grid and time array below.
CFL = 0.20;
dt = CFL * min(dx,dy) / c0;
tEnd = 18e-6;        % total simulation time [s]
Nt = round(tEnd/dt) + 1;
kgrid.setTime(Nt, dt);
t_array = kgrid.t_array;

%% 2) Define 16-element linear array at shallow depth
nElements = 16;
elementWidth = 0.40e-3; % element width [m]
kerf = 0.10e-3;         % gap between elements [m]
pitch = elementWidth + kerf;
elementWidthPts = max(1,round(elementWidth/dx));
pitchPts = max(elementWidthPts,round(pitch/dx));
ySource = 8;             % source row, away from boundary
centerX = round(Nx/2);

source.p_mask = false(Nx,Ny);
elementMasks = cell(nElements,1);
elementCenterIdx = zeros(nElements,1);

firstCenter = centerX - round(((nElements-1)*pitch)/(2*dx));
for e = 1:nElements
    xc = firstCenter + round((e-1)*pitch/dx);
    elementCenterIdx(e) = xc;
    x1 = max(1, xc-floor(elementWidthPts/2));
    x2 = min(Nx, x1+elementWidthPts-1);
    mask = false(Nx,Ny);
    mask(x1:x2,ySource) = true;
    elementMasks{e} = mask;
    source.p_mask = source.p_mask | mask;
end

% A short tone burst with a Hann envelope.
tone = sin(2*pi*f0*t_array);
nBurst = min(numel(t_array),round(nCycles/(f0*dt)));
window = zeros(size(t_array));
if nBurst >= 2
    window(1:nBurst) = 0.5*(1-cos(2*pi*(0:nBurst-1)/(nBurst-1)));
end
signal = sourcePressure * tone .* window;

% k-Wave orders source.p rows according to find(source.p_mask).
% Build one time series per active grid point, with optional element delays.
activeIdx = find(source.p_mask);
source.p = zeros(numel(activeIdx), Nt);
source.p_mode = 'additive';

% Steering angle: 0 = straight ahead. Positive angle steers laterally.
steeringAngleDeg = 0;
steeringAngle = deg2rad(steeringAngleDeg);
elementX_m = (elementCenterIdx - centerX) * dx;
elementDelay = elementX_m * sin(steeringAngle) / c0;
% Delay each element by shifting its burst in time.
for e = 1:nElements
    elementLinearIdx = find(elementMasks{e});
    [~, rowsInSource] = ismember(elementLinearIdx, activeIdx);
    rowsInSource = rowsInSource(rowsInSource > 0);
    delaySamples = round(elementDelay(e)/dt);
    for r = rowsInSource(:)'
        shifted = zeros(1,Nt);
        srcStart = max(1,1+delaySamples);
        srcEnd = min(Nt,nBurst+delaySamples);
        if srcStart <= srcEnd
            origStart = srcStart-delaySamples;
            origEnd = srcEnd-delaySamples;
            shifted(srcStart:srcEnd) = signal(origStart:origEnd);
        end
        source.p(r,:) = shifted;
    end
end

%% 3) Record acoustic pressure field
sensor.mask = true(Nx,Ny);
sensor.record = {'p_max','p_min','p_rms'};
input_args = {'PMLSize', 20, 'PlotSim', false, 'DataCast', 'single'};

% Run simulation. On some k-Wave versions, GPU-specific options differ.
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% Reshape recorded arrays back to Nx-by-Ny images.
pmax = reshape(sensor_data.p_max, Nx, Ny);
pmin = reshape(sensor_data.p_min, Nx, Ny);
prms = reshape(sensor_data.p_rms, Nx, Ny);

%% 4) Visualize pressure field and beam profiles
x_mm = kgrid.x_vec * 1e3;
y_mm = kgrid.y_vec * 1e3;

figure('Name','Maximum positive pressure');
imagesc(y_mm, x_mm, pmax/1e3); axis image; set(gca,'YDir','normal');
xlabel('Depth y (mm)'); ylabel('Lateral x (mm)');
title('Maximum positive pressure (kPa)'); colorbar;

figure('Name','Peak negative pressure');
imagesc(y_mm, x_mm, pmin/1e3); axis image; set(gca,'YDir','normal');
xlabel('Depth y (mm)'); ylabel('Lateral x (mm)');
title('Minimum pressure (kPa; more negative = stronger rarefaction)'); colorbar;

figure('Name','RMS pressure');
imagesc(y_mm, x_mm, prms/1e3); axis image; set(gca,'YDir','normal');
xlabel('Depth y (mm)'); ylabel('Lateral x (mm)');
title('RMS pressure (kPa)'); colorbar;

% Beam profile at a selected depth.
depth_mm = 8;
[~, depthIdx] = min(abs(y_mm-depth_mm));
figure('Name','Lateral beam profile');
plot(x_mm, pmax(:,depthIdx)/1e3, 'LineWidth',1.5);
xlabel('Lateral x (mm)'); ylabel('Peak positive pressure (kPa)');
title(sprintf('Lateral profile at depth %.2f mm', y_mm(depthIdx))); grid on;

% Basic summary metrics
[maxP, maxLinearIdx] = max(pmax(:));
[minP, minLinearIdx] = min(pmin(:));
[ixMax, iyMax] = ind2sub(size(pmax), maxLinearIdx);
[ixMin, iyMin] = ind2sub(size(pmin), minLinearIdx);

fprintf('Array: %d elements, f0 = %.2f MHz\n', nElements, f0/1e6);
fprintf('Nominal source pressure = %.1f kPa (not automatically calibrated to face pressure)\n', sourcePressure/1e3);
fprintf('Max positive pressure in domain = %.1f kPa at x=%.2f mm, y=%.2f mm\n', ...
    maxP/1e3, x_mm(ixMax), y_mm(iyMax));
fprintf('Min pressure in domain = %.1f kPa at x=%.2f mm, y=%.2f mm\n', ...
    minP/1e3, x_mm(ixMin), y_mm(iyMin));

%% 5) Suggested parameter studies
% Re-run while changing one parameter at a time:
% - f0: frequency (e.g., 0.5, 1, 2 MHz)
% - sourcePressure: calibrate using hydrophone measurement
% - steeringAngleDeg: beam steering
% - elementWidth, kerf, pitch: array geometry
% - per-element phase delays: focusing / apodization
% - medium.sound_speed and medium.density: medium properties
%
% To focus at (xf, yf), calculate each element's propagation distance:
%   d_e = sqrt((xf - x_e)^2 + yf^2)
%   delay_e = (max(d_e) - d_e)/c0
% Then shift each element's source signal by delay_e/dt samples.
%
% Note: the 2D model represents a line array in a 2D medium, not the full
% 3D field of a finite-elevation transducer. Use kspaceFirstOrder3D for
% elevation effects, side lobes, and realistic 3D beam shape.
