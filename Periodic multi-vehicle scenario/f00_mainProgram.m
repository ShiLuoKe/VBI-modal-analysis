%% f00_mainProgram.m
% =========================================================================
% Description : Main program for Periodic multi-vehicle scenario
% Author      : Shi LuoKe (史罗克)
% Date        : 2026-04-18
% -------------------------------------------------------------------------
% Acknowledgment: To improve code readability, this code was organized
%                 and formatted with the assistance of Claude (Anthropic).
%                 The correctness of the code has been verified by the author.
%                 
% -------------------------------------------------------------------------
% Note        : This program extends the single-vehicle VBI code to a
%               periodic multi-vehicle scenario. An idealized infinite
%               sequence of equally spaced vehicles travels across the
%               bridge at a constant speed, enabling eigenvalue analysis
%               of VBI systems under sustained vehicle loading.
%
%               Reference: XXXX (to be added after publication)
% =========================================================================

clc
close all
clear

%% Step 1: Load input parameters
sysParams = f01_loadInputParams();

%% Step 2: Build element matrices and compute uncoupled frequencies
[beamElem, theorFreq, sysParams] = f02_buildElementMatrices(sysParams);

%% Step 3: Assemble time-varying global matrices (multi-vehicle VBI)
globalMat = f03_assembleGlobalMatrices(sysParams, beamElem);

%% Step 4: Solve coupled eigenvalue problem at each time step
eigenResult = f04_solveEigenValue(globalMat, sysParams, theorFreq);

%% Step 5: Solve dynamic response using Newmark-beta method
dynResult = f05_solveNewmarkBeta(globalMat, sysParams);

%% Step 6: Post-processing
f06_postProcess(eigenResult, dynResult, theorFreq, sysParams);