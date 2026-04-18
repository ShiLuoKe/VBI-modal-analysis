%% f00_mainProgram.m
% =========================================================================
% Description : Main program for Single-vehicle scenario (gamma varies)
% Author      : Shi LuoKe (史罗克)
% Date        : 2026-04-08
% -------------------------------------------------------------------------
% Acknowledgment: To improve code readability, this code was organized
%                 and formatted with the assistance of Claude (Anthropic). 
%                 The correctness of the code has been verified by the author.
% -------------------------------------------------------------------------
% Note        : This program computes the dynamic response and eigenvalues
%               of a single-degree-of-freedom sprung mass moving on a
%               Bernoulli-Euler beam.
%
%               The primary objective is to demonstrate how different
%               configurations (controlled by gamma)
%               lead to variations in the system eigenvalues.
%
%               To modify gamma, please do so in the subfunction f01.
%
%               Reference: XXXX (to be added after publication)
% =========================================================================

clc             % Clear command window
close all       % Close all figures
clear           % Clear workspace

%% ================================================================
%  gamma = 1
%  Proposed value for modal anaylsis
%  To change gamma, modify the value in subfunction f01.
%% ================================================================

% Step 1: Load input parameters
sysParams = f01_loadInputParams();

% Step 2: Build element matrices and compute uncoupled frequencies
[beamElem, theorFreq, sysParams] = f02_buildElementMatrices(sysParams);

% Step 3: Assemble time-varying global matrices (including VBI element)
globalMat = f03_assembleGlobalMatrices(sysParams, beamElem);

% Step 4: Solve coupled eigenvalue problem at each time step
eigenResult = f04_solveEigenValue(globalMat, sysParams, theorFreq);

% Step 5: Solve dynamic response using Newmark-beta method
dynResult = f05_solveNewmarkBeta(globalMat, sysParams);
%%
% Step 6: Post-processing
f06_postProcess(eigenResult, dynResult, theorFreq, sysParams);