function [opt, pslDataNameOutput] = RunMission_Via_LineVis_Call(userInterface)	
	%% **********Note**********
	%% This Script is Dedicated to the Scenario that the Render "LineVis" (https://github.com/chrismile/LineVis)
	%%	is Used as the Frontend to Call the "3D-TSV", where the PSL Data is Written into "NAME_psl.dat" for the Info. Exchange
	global tracingFuncHandle_;
	global tracingStepWidth_;
	global majorPSLindexList_;
	global mediumPSLindexList_;
	global minorPSLindexList_;
	global nodeWiseStressField_;
	
	%%1. Initialize Experiment Environment
	%%1.1 Variable Declaration	
	tStart = tic;
	opt = 0; 
	pslDataNameOutput = [];
	GlobalVariables;	
	%%1.2 Decode Input Arguments
	fileName = userInterface.fileName;
	if ~strcmp(dataName_, fileName)		
		disp('Loading Dataset....');
		ImportStressFields(fileName);
		dataName_ = fileName;
	end
	
	lineDensCtrl = userInterface.lineDensCtrl;
	minimumEpsilon_ = min(boundingBox_(2,:)-boundingBox_(1,:))/lineDensCtrl;
	numLevels = userInterface.numLevels;
	seedStrategy = userInterface.seedStrategy;
	seedDensCtrl = userInterface.seedDensCtrl;
	
	selectedPrincipalStressField_ = userInterface.selectedPrincipalStressField;
	numPSF = length(selectedPrincipalStressField_);
	for ii=1:numPSF
		iPSF = selectedPrincipalStressField_(ii);
		switch iPSF
			case 1, PSLsAppearanceOrder_(end+1,:) = [1 0];
			case 2, PSLsAppearanceOrder_(end+1,:) = [2 0];
			case 3, PSLsAppearanceOrder_(end+1,:) = [3 0];
		end
	end
	
	mergingOpt_ = userInterface.mergingOpt;
	if ~mergingOpt_, numLevels = 1; end
	multiMergingThresholdsCtrl_ = userInterface.multiMergingThresholds;
	
	snappingOpt_ = userInterface.snappingOpt;
	
	permittedMaxAdjacentTangentAngleDeviation_ = userInterface.maxAngleDevi;
	
	traceAlg_ = userInterface.traceAlgorithm;
	if strcmp(meshType_, 'CARTESIAN_GRID')
		switch traceAlg_
			case 'Euler', tracingFuncHandle_ = @TracingPSL_Euler_CartesianMesh;
			case 'RK2', tracingFuncHandle_ = @TracingPSL_RK2_CartesianMesh;
			case 'RK4', tracingFuncHandle_ = @TracingPSL_RK4_CartesianMesh;
		end
		if ~nodeWiseStressField_
			tracingFuncHandle_ = @TracingPSL_Euler_CartesianMesh;
		end
	else
		switch traceAlg_
			case 'Euler', tracingFuncHandle_ = @TracingPSL_Euler_UnstructuredMesh;
			case 'RK2', tracingFuncHandle_ = @TracingPSL_RK2_UnstructuredMesh;
			case 'RK4', tracingFuncHandle_ = @TracingPSL_RK4_UnstructuredMesh;
		end
		if ~nodeWiseStressField_
			tracingFuncHandle_ = @TracingPSL_Euler_UnstructuredMesh;
		end
	end	

	%%2. Seeding
	GenerateSeedPoints(seedStrategy, seedDensCtrl);
	
	%%3. PSL Generation
	if strcmp(meshType_, 'CARTESIAN_GRID')
		tracingStepWidth_ = integratingStepScalingFac_ * eleSize_;
		integrationStepLimit_ = ceil(1.5*norm(boundingBox_(2,:)-boundingBox_(1,:))/tracingStepWidth_);
	else	
		tracingStepWidth_ = integratingStepScalingFac_ * eleSizeList_;
		integrationStepLimit_ = ceil(1.5*norm(boundingBox_(2,:)-boundingBox_(1,:))/median(tracingStepWidth_));
	end
	
	majorPSLindexList_ = struct('arr', []); majorPSLindexList_ = repmat(majorPSLindexList_, 1, numLevels);
	mediumPSLindexList_ = struct('arr', []); mediumPSLindexList_ = repmat(mediumPSLindexList_, 1, numLevels);     
	minorPSLindexList_ = struct('arr', []); minorPSLindexList_ = repmat(minorPSLindexList_, 1, numLevels);	
	index = 1;
	while index<=numLevels
		iEpsilon = minimumEpsilon_ * 2^(numLevels-index);
		GenerateSpaceFillingPSLs(iEpsilon);
		majorPSLindexList_(index).arr = 1:length(majorPSLpool_);
		mediumPSLindexList_(index).arr = 1:length(mediumPSLpool_);
		minorPSLindexList_(index).arr = 1:length(minorPSLpool_);
		index = index + 1;
	end	
	
    %%4. Building Hierarchy
    BuildPSLs4Hierarchy();
	
	%%5. Write&Print Results
	[~,~,fileExtension] = fileparts(dataName_);
	pslDataNameOutput = strcat(erase(dataName_, fileExtension), '_psl.dat');
	ExportResult(pslDataNameOutput);
	opt = 1;
	tEnd = toc(tStart);
	PrintAlgorithmStatistics(tEnd);
end