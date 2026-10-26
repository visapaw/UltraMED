# Ultrasound–Microbubble Simulation

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

- Rayleigh–Plesset (RP): อธิบายการเคลื่อนที่ของผิวฟองผ่านรัศมี \(R(t)\) โดยคำนึงถึงความดัน ความเฉื่อย ความหนืดของของเหลว และแรงตึงผิว
- Marmottant shell model: เพิ่มพฤติกรรมของ lipid shell เช่น การยุบตัว การยืดหยุ่น และสภาวะที่ shell ถือว่าแตก
- <img width="1344" height="286" alt="image" src="https://github.com/user-attachments/assets/72de2646-7944-4b42-802e-1448c6e7f166" />
- <img width="1324" height="262" alt="image" src="https://github.com/user-attachments/assets/3720f73b-2b07-43a7-b2a4-5b13f07ff9b6" />

### 2. Transducer Acoustic Field


## Files

| File | Description |
|---|---|
| `bubble_dynamics_marmottant.m` | Bubble oscillation and modeled shell-rupture analysis |
| `kwave_16element_linear_array.m` | 2D acoustic-field simulation for a 16-element array |


## Limitations
ปัญหาคือยังไม่มีค่าk (surface dilatational viscosity), X (compression modulus), Sigma0 ที่เฉพาะเจาะจงสำหรับการใช้lipid shell ratio()และAMC
