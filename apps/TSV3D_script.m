%% 3D-TSV
%%The 3D Trajectory-based Stress Visualizer (3D-TSV), a visual analysis tool for the exploration 
%%of the principal stress directions in 3D solids under load.

%%This repository was created for the paper "3D-TSV: The 3D Trajectory-based Stress Visualizer" 
%%	by Junpeng Wang, Christoph Neuhauser, Jun Wu, Xifeng Gao and Rüdiger Westermann, 
%%which was submitted to the journal "Advances in Engineering Software", and also available in arXiv (2112.09202)

%% ============= Test data sets can be found in =============
%%	https://syncandshare.lrz.de/getlink/fi4W4EGjZSzMzCvxkEf9L3Aw/

clc;
addpath('../backend');
userInterface = InterfaceStruct();
%% Prepare the data set, and uncomment one of the experiments below to run the 3D-TSV for testing the examples presented in the paper
%% *****please be sure to relate the correct directory of data set*****

%% Some Examples used in the paper
%% =======================================kitten=======================================
userInterface.fileName = '../data/kitten.TSV';
userInterface.lineDensCtrl = 20;
userInterface.numLevels = 1;
userInterface.seedStrategy = 'Volume';
userInterface.seedDensCtrl = 2;
userInterface.selectedPrincipalStressField = [3];
userInterface.mergingOpt = 1;
userInterface.snappingOpt = 0;
userInterface.maxAngleDevi = 6;
userInterface.traceAlgorithm = 'RK2';


%% =======================================cantilever=======================================
% userInterface.fileName = '../data/cantilever3D.TSVcarti';
% userInterface.lineDensCtrl = 5; 
% userInterface.numLevels = 1;
% userInterface.seedStrategy = 'Volume';
% userInterface.seedDensCtrl = 5; 
% userInterface.selectedPrincipalStressField = [1, 3];
% userInterface.mergingOpt = 1;
% userInterface.snappingOpt = 0;
% userInterface.maxAngleDevi = 20;
% userInterface.traceAlgorithm = 'Euler';


%% =======================================arched_bridge=======================================
%% ---------------------------------------Experiment 1: (fig. 3)
% userInterface.fileName = '../data/arched_bridge.TSVcarti';
% userInterface.lineDensCtrl = 6; 
% userInterface.numLevels = 1;
% userInterface.seedStrategy = 'Volume';
% userInterface.seedDensCtrl = 2; 
% userInterface.selectedPrincipalStressField = [1, 2, 3];
% userInterface.mergingOpt = 1;
% userInterface.snappingOpt = 0;
% userInterface.maxAngleDevi = 20; 
% userInterface.traceAlgorithm = 'RK2';
%% ---------------------------------------Experiment 2: (fig.4)
% userInterface.fileName = '../data/arched_bridge.TSVcarti';
% userInterface.lineDensCtrl = 24; 
% userInterface.numLevels = 4;
% userInterface.seedStrategy = 'Volume';
% userInterface.seedDensCtrl = 2; 
% userInterface.selectedPrincipalStressField = [1, 3];
% userInterface.mergingOpt = 1;
% userInterface.snappingOpt = 0;
% userInterface.maxAngleDevi = 20; 
% userInterface.traceAlgorithm = 'RK2';


%% =======================================rod=======================================
% userInterface.fileName = '../data/rod.carti';
% userInterface.lineDensCtrl = 5;
% userInterface.numLevels = 1;
% userInterface.seedStrategy = 'Volume';
% userInterface.seedDensCtrl = 3;
% userInterface.selectedPrincipalStressField = [1, 3];
% userInterface.mergingOpt = 1;
% userInterface.snappingOpt = 0;
% userInterface.maxAngleDevi = 20;
% userInterface.traceAlgorithm = 'RK2';


%% =======================================femur=======================================
% userInterface.fileName = '../data/femur.TSVcarti';
% userInterface.lineDensCtrl = 18;
% userInterface.numLevels = 3;
% userInterface.seedStrategy = 'Volume';
% userInterface.seedDensCtrl = 4;
% userInterface.selectedPrincipalStressField = [1, 2, 3];
% userInterface.mergingOpt = 1;
% userInterface.snappingOpt = 0;
% userInterface.maxAngleDevi = 10;
% userInterface.traceAlgorithm = 'RK2';


%% =======================================bracket=======================================
% userInterface.fileName = '../data/bracket.TSVcarti';
% userInterface.lineDensCtrl = 12;
% userInterface.numLevels = 3;
% userInterface.seedStrategy = 'Volume';
% userInterface.seedDensCtrl = 4;
% userInterface.selectedPrincipalStressField = [1, 3];
% userInterface.mergingOpt = 1;
% userInterface.snappingOpt = 0;
% userInterface.maxAngleDevi = 6;
% userInterface.traceAlgorithm = 'RK2';


%% =======================================parts1=======================================
% userInterface.fileName = '../data/parts1.TSV';
% userInterface.lineDensCtrl = 30;
% userInterface.numLevels = 3;
% userInterface.seedStrategy = 'Volume';
% userInterface.seedDensCtrl = 5;
% userInterface.selectedPrincipalStressField = [1, 3];
% userInterface.mergingOpt = 1;
% userInterface.snappingOpt = 0;
% userInterface.maxAngleDevi = 6;
% userInterface.traceAlgorithm = 'RK2';

%% =======================================bearing=======================================
% userInterface.fileName = '../data/bearing.TSV';
% userInterface.lineDensCtrl = 21;
% userInterface.numLevels = 3;
% userInterface.seedStrategy = 'Volume';
% userInterface.seedDensCtrl = 4;
% userInterface.selectedPrincipalStressField = [1, 3];
% userInterface.mergingOpt = 1;
% userInterface.snappingOpt = 0;
% userInterface.maxAngleDevi = 6;
% userInterface.traceAlgorithm = 'RK2';

RunMission(userInterface);

%%PSLs Visualization
%% ======Syntax======
% DrawPSLs(imOpt, imVal, pslGeo, stressComponentOpt, lw, minLength);
DrawPSLs(["Geo", "Geo", "Geo"], [0,0,0], ["TUBE", "TUBE", "TUBE"], 'None', 0.5, 20);
%% Show if Necessary
% DrawSeedPoints(0.5, 'inputSeeds');
