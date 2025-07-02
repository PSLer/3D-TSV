%% This script declares some global variables for this program

%% 1. Dataset info.
global dataName_; 
global meshType_;
global boundingBox_; 
global nelx_; 
global nely_; 
global nelz_; 
global numNodes_;
global nodeCoords_; 
global numEles_; 
global eNodMat_;
global nodState_; 
global nodStruct_; %% Adjacent Elements for each Mesh Vertex (only for Unstructured Hex-mesh)
global eleStruct_; %% Element, face normals, adjacent element of each face 
global cartesianStressField_; %% at each row in the oder: Sigma_xx, Sigma_yy, Sigma_zz, Sigma_yz, Sigma_zx, Sigma_xy
global loadingCond_; 
global fixingCond_; 	
global eleSize_; %% Element Size for Cartesian Mesh
global eleSizeList_; %% Element Sizes for Unstructured Hex-mesh
global integrationStepLimit_; %% Maximum Integration Steps
global silhouetteStruct_; %% Patches for Draw Silhoutte of Stress Field
global carNodMapBack_; %% Only for Cartesian Mesh
global surfaceMeshNodeCoords_;
global surfaceMeshElements_;

%% 2. Algorithm Control
%% 2.0 Selected Principal Stress Fields to be Visualized
global selectedPrincipalStressField_; 
%% 2.1 Integrating Step Scaling Factor, positively related to the integrating step size
global integratingStepScalingFac_; integratingStepScalingFac_ = 1;
%% 2.2 %% Tracing PSL stops when the angle deviation between the neighboring tangents is larger than 180/permittedMaxAdjacentTangentAngleDeviation_
global permittedMaxAdjacentTangentAngleDeviation_;
%% 2.3 Store the Original Seed Points
global seedPointsHistory_;
%% 2.4 Control Seed Point Density under 'Volume' Seeding Strategy
global seedSpan4VolumeOptCartesianMesh_; 
%% 2.5 Relaxing Merge Operation via Epsilon*relaxedFactor_
global relaxedFactor_; relaxedFactor_ = 1;
global multiMergingThresholdsCtrl_;
%% 2.6 Snapping PSLs (1) or not (0) when they are too close
global snappingOpt_; 
%% 2.7 Excluding PSLs with less than 'minLengthVisiblePSLs_' Integrating Steps
global minLengthVisiblePSLs_; minLengthVisiblePSLs_ = 5;
%% 2.8 Merging or not. FALSE == Generating PSLs Brutally from Selected Seed Points 
global mergingOpt_; 
%% 2.9 PSL Tracing Algorithm
global traceAlg_; 
%% 2.10 Minimum Merging Threshold 
global minimumEpsilon_;

%% 3. Result
global majorPSLpool_; majorPSLpool_ = PrincipalStressLineStruct();
global mediumPSLpool_; mediumPSLpool_ = PrincipalStressLineStruct();
global minorPSLpool_; minorPSLpool_ = PrincipalStressLineStruct();
global majorCoordList_; majorCoordList_ = [];
global mediumCoordList_; mediumCoordList_ = [];
global minorCoordList_; minorCoordList_ = [];
global majorHierarchy_; majorHierarchy_ = [];
global mediumHierarchy_; mediumHierarchy_ = [];
global minorHierarchy_; minorHierarchy_ = [];
global PSLsAppearanceOrder_; PSLsAppearanceOrder_ = [];

%% 4. Visualization (only for Matlab Demo Code)
global lineWidth_; lineWidth_ = 1;
global silhouetteOpacity_; silhouetteOpacity_ = 0.1;
global tubeShapedPSLs_Patches_; tubeShapedPSLs_Patches_ = [];
global ribbonShapedPSLs_Patches_; ribbonShapedPSLs_Patches_ = [];