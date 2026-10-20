# Ultrasound–Microbubble Simulation
This project investigates the effects of ultrasound parameters on lipid-shell microbubble dynamics and explores acoustic-field simulation for a 16-element soft-piezoceramic transducer.

## Project Overview

The microbubble system consists of:
- Lipid shell
- Fluorescent surrogate: AMC (7-Amino-4-methylcoumarin)
- Gas core: Nitrogen (N₂)

Transducer:
16-element soft-piezoceramic

## Goals

1. Ultrasound ที่ความถี่ไหนทำให้ microbubble แตก
2. ออกแบบ 16-element soft-piezoceramic ที่ให้ค่าที่เหมาะสม

## Equation used and Referrences

### 1. Bubble Dynamics


### 2. Transducer Acoustic Field


## Files

| File | Description |
|---|---|
| `bubble_dynamics_marmottant.m` | Bubble oscillation and modeled shell-rupture analysis |
| `kwave_16element_linear_array.m` | 2D acoustic-field simulation for a 16-element array |


## Limitations

The simulation parameters are initial illustrative values, not validated measurements of the actual formulation. The effects of AMC on shell properties must be established experimentally. The transducer model does not directly simulate piezoelectric coupling.
