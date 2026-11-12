% Copyright (c) 2025-2026 Dennis Huber, Celina Bräutigam
% SPDX-License-Identifier: LicenseRef-PolyForm-Noncommercial-1.0.0
%% Setup
% Get the Test directory
test_dir = fileparts(mfilename('fullpath'));

% Change to that directory
cd(test_dir);

% Create output directories if they do not exist
if ~exist("TestBipartiteComplete", 'dir')
    mkdir("TestBipartiteComplete");
end
if ~exist("TestTripartiteComplete", 'dir')
    mkdir("TestTripartiteComplete");
end
if ~exist("TestStackedPlot", 'dir')
    mkdir("TestStackedPlot");
end


%% Bipartite Tests: Werner stae
rho = 0.3/4 * eye(4) + 0.7 * 1/2*[1 0 0 1]'*[1 0 0 1];

% Complete View
path = "TestBipartiteComplete/Werner.pdf";
visualize(rho, path, "BeadTypes", "M", "CameraAngle",[47, 30], ...
    "GridPoints",[200 400], "Id", false, "AutoScale", true, "EntanglementArcs", true, "FontSize",12);

% All possible bead types tested in one run
path = "TestBipartiteComplete/WernerAll.pdf";
visualize(rho, path, "BeadTypes", ["M", "T", "E", "C", "Sigma", "Sigma_Cl", "Sigma_D", "D", "M_Cl"], "CameraAngle",[47, 30], ...
    "GridPoints",[200 400], "Id", false, "AutoScale", true, "EntanglementArcs", true, "FontSize",12, "SingleView", true);

% All bead types tested individually
% T-Beads
path = "TestBipartiteComplete/WernerT.pdf";
visualize(rho, path, "BeadTypes", "T", "CameraAngle",[47, 30], ...
    "GridPoints",[200 400], "Id", false, "AutoScale", true, "EntanglementArcs", true, "FontSize",12, "SingleView", true);

% E-Beads
path = "TestBipartiteComplete/WernerE.pdf";
visualize(rho, path, "BeadTypes", "E", "CameraAngle",[47, 30], ...
    "GridPoints",[200 400], "Id", false, "AutoScale", true, "EntanglementArcs", true, "FontSize",12, "SingleView", true);

% C-Beads
path = "TestBipartiteComplete/WernerC.pdf";
visualize(rho, path, "BeadTypes", "C", "CameraAngle",[47, 30], ...
    "GridPoints",[200 400], "Id", false, "AutoScale", true, "EntanglementArcs", true, "FontSize",12, "SingleView", true);

% M-Beads
path = "TestBipartiteComplete/WernerM.pdf";
visualize(rho, path, "BeadTypes", "M", "CameraAngle",[47, 30], ...
    "GridPoints",[200 400], "Id", false, "AutoScale", true, "EntanglementArcs", true, "FontSize",12, "SingleView", true);

% Sigma-Beads
path = "TestBipartiteComplete/WernerSigma.pdf";
visualize(rho, path, "BeadTypes", "Sigma", "CameraAngle",[47, 30], ...
    "GridPoints",[200 400], "Id", false, "AutoScale", true, "EntanglementArcs", true, "FontSize",12, "SingleView", true);

% M_Cl-Beads
path = "TestBipartiteComplete/WernerM_Cl.pdf";
visualize(rho, path, "BeadTypes", "M_Cl", "CameraAngle",[47, 30], ...
    "GridPoints",[200 400], "Id", false, "AutoScale", true, "EntanglementArcs", true, "FontSize",12, "SingleView", true);

% D-Beads
path = "TestBipartiteComplete/WernerD.pdf";
visualize(rho, path, "BeadTypes", "D", "CameraAngle",[47, 30], ...
    "GridPoints",[200 400], "Id", false, "AutoScale", true, "EntanglementArcs", true, "FontSize",12, "SingleView", true);

% Sigma_Cl-Beads
path = "TestBipartiteComplete/WernerSigma_Cl.pdf";
visualize(rho, path, "BeadTypes", "Sigma_Cl", "CameraAngle",[47, 30], ...
    "GridPoints",[200 400], "Id", false, "AutoScale", true, "EntanglementArcs", true, "FontSize",12, "SingleView", true);

% Sigma_D-Beads
path = "TestBipartiteComplete/WernerSigma_D.pdf";
visualize(rho, path, "BeadTypes", "Sigma_D", "CameraAngle",[47, 30], ...
    "GridPoints",[200 400], "Id", false, "AutoScale", false, "EntanglementArcs", true, "FontSize",12, "SingleView", true);



%% GHZ Werner
rho = 0.3/8 * eye(8) + 0.7 * 1/2*[1 0 0 0 0 0 0 1]'*[1 0 0 0 0 0 0 1];

% All possible bead types tested in one run
path = "TestTripartiteComplete/GHZWerner.pdf";
visualize(rho, path, "BeadTypes", ["M", "T", "E", "C", "Sigma", "Sigma_Cl", "Sigma_D", "D", "M_Cl"], "CameraAngle",[54, 45], ...
    "GridPoints",[200 400],"TargetQubits",[1, 1, 2, 1, 2], "Id", false, "AutoScale", true, "EntanglementArcs", true, "PadVectors", 0, "FontSize",12);

% All bead types tested individually
% T-Beads
path = "TestTripartiteComplete/GHZWernerT_.pdf";
visualize(rho, path, "BeadTypes", "T", "CameraAngle",[54, 45], ...
    "GridPoints",[200 400],"TargetQubits",[1, 1, 2, 1, 2], "Id", false, "AutoScale", true, "EntanglementArcs", true, "PadVectors", 0, "FontSize",12);

% E-Beads
path = "TestTripartiteComplete/GHZWernerE_.pdf";
visualize(rho, path, "BeadTypes", "E", "CameraAngle",[54, 45], ...
    "GridPoints",[200 400],"TargetQubits",[1, 1, 2, 1, 2], "Id", false, "AutoScale", true, "EntanglementArcs", true, "PadVectors", 0, "FontSize",12);

% C-Beads
path = "TestTripartiteComplete/GHZWernerC_.pdf";
visualize(rho, path, "BeadTypes", "C", "CameraAngle",[54, 45], ...
    "GridPoints",[200 400],"TargetQubits",[1, 1, 2, 1, 2], "Id", false, "AutoScale", true, "EntanglementArcs", true, "PadVectors", 0,"FontSize",12);

% M-Beads
path = "TestTripartiteComplete/GHZWernerM_.pdf";
visualize(rho, path, "BeadTypes", "M", "CameraAngle",[54, 45], ...
    "GridPoints",[200 400],"TargetQubits",[1, 1, 2, 1, 2], "Id", false, "AutoScale", true, "EntanglementArcs", true, "PadVectors", 0, "FontSize",12);

% Sigma-Beads
path = "TestTripartiteComplete/GHZWernerSigma_.pdf";
visualize(rho, path, "BeadTypes", "Sigma", "CameraAngle",[54, 45], ...
    "GridPoints",[200 400],"TargetQubits",[1, 1, 2, 1, 2], "Id", false, "AutoScale", true, "EntanglementArcs", true, "PadVectors", 0, "FontSize",12);

% M_Cl-Beads
path = "TestTripartiteComplete/GHZWernerM_Cl_.pdf";
visualize(rho, path, "BeadTypes", "M_Cl", "CameraAngle",[54, 45], ...
    "GridPoints",[200 400],"TargetQubits",[1, 1, 2, 1, 2], "Id", false, "AutoScale", true, "EntanglementArcs", true, "PadVectors", 0, "FontSize",12);

% D-Beads
path = "TestTripartiteComplete/GHZWernerD_.pdf";
visualize(rho, path, "BeadTypes", "D", "CameraAngle",[54, 45], ...
    "GridPoints",[200 400],"TargetQubits",[1, 1, 2, 1, 2], "Id", false, "AutoScale", true, "EntanglementArcs", true, "PadVectors", 0, "FontSize",12);

% Sigma_Cl-Beads
path = "TestTripartiteComplete/GHZWernerSigma_Cl_.pdf";
visualize(rho, path, "BeadTypes", "Sigma_Cl", "CameraAngle",[54, 45], ...
    "GridPoints",[200 400],"TargetQubits",[1, 1, 2, 1, 2], "Id", false, "AutoScale", true, "EntanglementArcs", true, "PadVectors", 0, "FontSize",12);

% Sigma_D-Beads
path = "TestTripartiteComplete/GHZWernerSigma_D_.pdf";
visualize(rho, path, "BeadTypes", "Sigma_D", "CameraAngle",[54, 45], ...
    "GridPoints",[200 400],"TargetQubits",[1, 1, 2, 1, 2], "Id", false, "AutoScale", true, "EntanglementArcs", true, "PadVectors", 0, "FontSize",12);

%% Stacked plots test
% Zurek example classical correlation -> Bell state
State00 = [1;0;0;0];
State11 = [0;0;0;1];

% Direct function handle
state = @(z) 0.5*(State00*State00'+State11*State11')+...
      (z/2)*(State00*State11'+State11*State00');
path = "TestStackedPlot/Zurek.pdf";
stacked_plot(state, path, "ParamValues", 0:0.25:1);

% Parsed function handle
state = "k*dens(|Psi+>) + (1-k)/4*eye(4)";
path = "TestStackedPlot/Werner.pdf";
stacked_plot(state, path, "ParamValues", 0:0.25:1, "Parameter","k");

%% W-Test for optimized separable structure
% No optimization
path = "TestTripartiteComplete/W_Clnonopt.pdf";
visualize("|W>", path, "BeadTypes", ["M", "M_Cl"], "CameraAngle",[54, 45], ...
    "GridPoints",[200 400],"TargetQubits",[1, 1, 2, 1, 2], "Id", false, "AutoScale", true, "EntanglementArcs", true, "PadVectors", 0, "FontSize",12);

% With optimization
path = "TestTripartiteComplete/W_Clopt.pdf";
visualize("|W>", path, "BeadTypes", ["M", "M_Cl"], "CameraAngle",[54, 45], ...
    "GridPoints",[200 400],"TargetQubits",[1, 1, 2, 1, 2], "Id", false, "AutoScale", true, "EntanglementArcs", true, "PadVectors", 0, "FontSize",12, "OptStructure", 1);
