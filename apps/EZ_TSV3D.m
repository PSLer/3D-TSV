%% This code is created for performing the 3D-TSV to the stress tensor field simulated on the 1st-order Quadrilateral and Triangular Mesh 
% Easy-Use Version
% Author: Junpeng Wang (junpeng.wang@tum.de)
% Date: 2023-12-12
clear all; clc;

global majorPSLpool_;
global mediumPSLpool_; 
global minorPSLpool_; 


%%1. Import Data
stressfileName = '../data/demoData_3D_Tet_femur.stress';
ImportStressFields(stressfileName);
figure; ShowProblemDescription();

%%2. PSLs generation
PSLsDensityCtrl = 10; %%The larger this value is, the denser the PSLs are.
TSV3D(PSLsDensityCtrl,3);
figure; ShowPSLs(100);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TSV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TSV3D(PSLsDensityCtrl, varargin)
	global boundingBox_;
	global eleCentroidList_;
	global numEles_;
	global eleSizeList_;
	
	global mergingThreshold_;
	global tracingStepWidth_;
	global integrationStepLimit_;
	global permittedMaxAdjacentTangentAngleDeviation_;
	global relaxedFactor_;
	
	global seedPointsHistory_;
	global seedAssociatedEles_;
	global seedPoints_;
	global seedPointsValence_;	
	global majorPSLpool_;
	global mediumPSLpool_;
    global minorPSLpool_; 
	global majorCoordList_; 
    global minorCoordList_;
	
	%%Settings
	mergingThreshold_ = min(boundingBox_(2,:)-boundingBox_(1,:))/PSLsDensityCtrl;
	permittedMaxAdjacentTangentAngleDeviation_ = 10;
	tracingStepWidth_ = 1.0 * min(boundingBox_(2,:)-boundingBox_(1,:))/100;
	tracingStepWidth_ = 0.5 * eleSizeList_;
	integrationStepLimit_ = ceil(1.5*norm(boundingBox_(2,:)-boundingBox_(1,:))/median(tracingStepWidth_));
	relaxedFactor_ = 1.0;
	
    if 1==nargin
        seedDensityCtrl = 1;
    else
        seedDensityCtrl = varargin{1};
    end
	seedAssociatedEles_ = 1:seedDensityCtrl:numEles_; seedAssociatedEles_ = seedAssociatedEles_(:);
	% seedAssociatedEles_(boundaryElements_) = [];
	seedPointsHistory_ = eleCentroidList_(seedAssociatedEles_,:);
	seedPointsValence_ = ones(size(seedPointsHistory_));
	numSeedPoints = size(seedPointsHistory_,1);
	seedPoints_ = [seedAssociatedEles_ seedPointsHistory_];
	
	InitializeMergingThresholdMap();
	
	majorPSLpool_ = PrincipalStressLineStruct();
	mediumPSLpool_ = PrincipalStressLineStruct();
	minorPSLpool_ = PrincipalStressLineStruct();

	%% Exclude the irrelated Principal Stress Fields
	selectedPrincipalStressField_ = [1 2 3];
	numPSF = length(selectedPrincipalStressField_);
	for ii=1:numPSF
		iPSF = selectedPrincipalStressField_(ii);
		switch iPSF
			case 1, seedPointsValence_(:,1) = 0; 
			case 2, seedPointsValence_(:,2) = 0; 
			case 3, seedPointsValence_(:,3) = 0;
		end
	end

	PreprocessSeedPoints();
	
	majorCoordList_ = [];
	mediumCoordList_ = [];
	minorCoordList_ = [];
	
	startCoord_ = boundingBox_(1,:) + sum(boundingBox_, 1)/2;
	
	%% Seeding
	its = 0;
	looper = sum(sum(seedPointsValence_));
	while looper<3*numSeedPoints
		its = its + 1;
		valenceMetric = sum(seedPointsValence_,2);
		%% 1st Priority: semi-empty seeds > empty seeds, which helps get PSLs intersection
		%% 2nd Priority: seeds with same valence, the one closest to the start point goes first	
		switch numPSF
			case 1
				unFinishedSppsValence2 = find(2==valenceMetric);
				[~, tarPos] = min(vecnorm(startCoord_-seedPoints_(unFinishedSppsValence2,end-2:end),2,2));
				spp = unFinishedSppsValence2(tarPos);
			case 2
				unFinishedSppsValence2 = find(2==valenceMetric); 
				if ~isempty(unFinishedSppsValence2) %% 1st Priority
					[~, tarPos] = min(vecnorm(startCoord_-seedPoints_(unFinishedSppsValence2,end-2:end),2,2)); %% 2nd Priority
					spp = unFinishedSppsValence2(tarPos);
				else
					unFinishedSppsValence1 = find(1==valenceMetric);
					[~, tarPos] = min(vecnorm(startCoord_-seedPoints_(unFinishedSppsValence1,end-2:end),2,2)); %% 2nd Priority
					spp = unFinishedSppsValence1(tarPos);			
				end					
			case 3
				unFinishedSppsValence12 = find(3>valenceMetric); 
				unFinishedSppsValence12 = unFinishedSppsValence12(valenceMetric(unFinishedSppsValence12)>0); 
				if ~isempty(unFinishedSppsValence12) %% 1st Priority
					[~, tarPos] = min(vecnorm(startCoord_-seedPoints_(unFinishedSppsValence12,end-2:end),2,2)); %% 2nd Priority
					spp = unFinishedSppsValence12(tarPos);
				else
					unFinishedSppsValence0 = find(0==valenceMetric);
					[~, tarPos] = min(vecnorm(startCoord_-seedPoints_(unFinishedSppsValence0,end-2:end),2,2)); %% 2nd Priority
					spp = unFinishedSppsValence0(tarPos);		
				end						
		end
		valences = seedPointsValence_(spp,:);						
		seed = seedPoints_(spp,:);		
		
		if 0==valences(1)
			seedPointsValence_(spp,1) = 1;
			majorPSL = Have1morePSL(seed, 'MAJOR');		
			if 0==majorPSL.length
				looper = sum(sum(seedPointsValence_)); 
				disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
					' Total.: ' sprintf('%6i',3*numSeedPoints)]);
				continue; 
			end			
			majorPSLpool_(end+1,1) = majorPSL;				
			majorCoordList_(end+1:end+majorPSL.length,:) = majorPSL.phyCoordList;
			sppsEmptyMajorValence = find(0==seedPointsValence_(:,1));
			if ~isempty(sppsEmptyMajorValence)
				[potentialDisListMajor, potentialPosListMajor] = GetDisListOfPointList2Curve(seedPoints_(...
						sppsEmptyMajorValence,:), [majorPSL.eleIndexList majorPSL.phyCoordList], 'MAJOR');					
				potentialSolidSppsMajor = find(potentialDisListMajor<relaxedFactor_);
				if ~isempty(potentialSolidSppsMajor)
					spps2BeMerged = sppsEmptyMajorValence(potentialSolidSppsMajor);
					seedPoints_(spps2BeMerged,:) = potentialPosListMajor(potentialSolidSppsMajor,:);								
					seedPointsValence_(spps2BeMerged,1) = 1;
					modifiedMediumValences = HighCurvatureModification(spps2BeMerged, 'MEDIUM');
					seedPointsValence_(modifiedMediumValences,2) = 1;					
					modifiedMinorValences = HighCurvatureModification(spps2BeMerged, 'MINOR');				
					seedPointsValence_(modifiedMinorValences,3) = 1;	
				end
			end				
		end		

		if 0==valences(2)
			seedPointsValence_(spp,2) = 1;
			mediumPSL = Have1morePSL(seed, 'MEDIUM');		
			if 0==mediumPSL.length
				looper = sum(sum(seedPointsValence_)); 
				disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
					' Total.: ' sprintf('%6i',3*numSeedPoints)]);
				continue; 
			end
			mediumPSLpool_(end+1,1) = mediumPSL;
			mediumCoordList_(end+1:end+mediumPSL.length,:) = mediumPSL.phyCoordList;
			sppsEmptyMediumValence = find(0==seedPointsValence_(:,2));
			if ~isempty(sppsEmptyMediumValence)
				[potentialDisListMedium, potentialPosListMedium] = GetDisListOfPointList2Curve(seedPoints_(...
						sppsEmptyMediumValence,:), [mediumPSL.eleIndexList mediumPSL.phyCoordList], 'MEDIUM');					
				potentialSolidSppsMedium = find(potentialDisListMedium<relaxedFactor_);
				if ~isempty(potentialSolidSppsMedium)
					spps2BeMerged = sppsEmptyMediumValence(potentialSolidSppsMedium);
					seedPoints_(spps2BeMerged,:) = potentialPosListMedium(potentialSolidSppsMedium,:);								
					seedPointsValence_(spps2BeMerged,2) = 1;
					modifiedMajorValences = HighCurvatureModification(spps2BeMerged, 'MAJOR');
					seedPointsValence_(modifiedMajorValences,1) = 1;
					modifiedMinorValences = HighCurvatureModification(spps2BeMerged, 'MINOR');
					seedPointsValence_(modifiedMinorValences,3) = 1;
				end
			end				
		end

		if 0==valences(3)
			seedPointsValence_(spp,3) = 1;			
			minorPSL = Have1morePSL(seed, 'MINOR');			
			if 0==minorPSL.length
				looper = sum(sum(seedPointsValence_)); 
				disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
					' Total.: ' sprintf('%6i',3*numSeedPoints)]);
				continue; 
			end		
			minorPSLpool_(end+1,1) = minorPSL;
			minorCoordList_(end+1:end+minorPSL.length,:) = minorPSL.phyCoordList;		
			sppsEmptyMinorValence = find(0==seedPointsValence_(:,3));
			if ~isempty(sppsEmptyMinorValence)   
				[potentialDisListMinor, potentialPosListMinor] = GetDisListOfPointList2Curve(seedPoints_(...
						sppsEmptyMinorValence,:), [minorPSL.eleIndexList minorPSL.phyCoordList], 'MINOR');					
				potentialSolidSppsMinor = find(potentialDisListMinor<relaxedFactor_);
				if ~isempty(potentialSolidSppsMinor)
					spps2BeMerged = sppsEmptyMinorValence(potentialSolidSppsMinor);
					seedPoints_(spps2BeMerged,:) = potentialPosListMinor(potentialSolidSppsMinor,:);
					seedPointsValence_(spps2BeMerged,3) = 1;				
					modifiedMajorValences = HighCurvatureModification(spps2BeMerged, 'MAJOR');					
					seedPointsValence_(modifiedMajorValences,1) = 1;
					modifiedMediumValences = HighCurvatureModification(spps2BeMerged, 'MEDIUM');
					seedPointsValence_(modifiedMediumValences,2) = 1;					
				end
			end					
		end		
		looper = sum(sum(seedPointsValence_));
		disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
			' Total.: ' sprintf('%6i',3*numSeedPoints)]);			
	end

	majorPSLpool_ = CompactStreamlines(majorPSLpool_, 5);
	mediumPSLpool_ = CompactStreamlines(mediumPSLpool_, 5);
	minorPSLpool_ = CompactStreamlines(minorPSLpool_, 5);	
end

function InitializeMergingThresholdMap()
	global numEles_;
	global mergingThreshold_;
	global mergingThresholdMap_;
	mergingThresholdMap_ = repmat(mergingThreshold_, numEles_, 3);
end

function PreprocessSeedPoints()
	global seedPoints_;
	global seedPointsValence_;
	global majorPSLpool_; 
	global mediumPSLpool_; 
	global minorPSLpool_; 
	global relaxedFactor_;
	global multiMergingThresholdsCtrl_;
	
	numMajorPSLs = length(majorPSLpool_);
	for ii=1:numMajorPSLs
		majorPSL = majorPSLpool_(ii);
		if majorPSL.length>0					
			sppsEmptyMajorValence = find(0==seedPointsValence_(:,1));
            if ~isempty(sppsEmptyMajorValence)
				[potentialDisListMajor, potentialPosListMajor] = GetDisListOfPointList2Curve(...	
					seedPoints_(sppsEmptyMajorValence,:), [majorPSL.eleIndexList majorPSL.phyCoordList], 'MAJOR');
				potentialSolidSppsMajor = find(potentialDisListMajor<=multiMergingThresholdsCtrl_(1)*relaxedFactor_);
				if ~isempty(potentialSolidSppsMajor)
					spps2BeMerged = sppsEmptyMajorValence(potentialSolidSppsMajor);							
					seedPoints_(spps2BeMerged,:) = potentialPosListMajor(potentialSolidSppsMajor,:);
					seedPointsValence_(spps2BeMerged,1) = 1;						
					modifiedMediumValences = HighCurvatureModification(spps2BeMerged, 'MEDIUM');
					seedPointsValence_(modifiedMediumValences,2) = 1;					
					modifiedMinorValences = HighCurvatureModification(spps2BeMerged, 'MINOR');
					seedPointsValence_(modifiedMinorValences,3) = 1;							
				end
			end
		end
	end

	numMediumPSLs = length(mediumPSLpool_);
	for ii=1:numMediumPSLs
		mediumPSL = mediumPSLpool_(ii);
		if mediumPSL.length>0					
			sppsEmptyMediumValence = find(0==seedPointsValence_(:,2));
            if ~isempty(sppsEmptyMediumValence)
				[potentialDisListMedium, potentialPosListMedium] = GetDisListOfPointList2Curve(...	
					seedPoints_(sppsEmptyMediumValence,:), [mediumPSL.eleIndexList mediumPSL.phyCoordList], 'MEDIUM');
				potentialSolidSppsMedium = find(potentialDisListMedium<=multiMergingThresholdsCtrl_(2)*relaxedFactor_);
				if ~isempty(potentialSolidSppsMedium)
					spps2BeMerged = sppsEmptyMediumValence(potentialSolidSppsMedium);							
					seedPoints_(spps2BeMerged,:) = potentialPosListMedium(potentialSolidSppsMedium,:);
					seedPointsValence_(spps2BeMerged,2) = 1;
					modifiedMajorValences = HighCurvatureModification(spps2BeMerged, 'MAJOR');
					seedPointsValence_(modifiedMajorValences,1) = 1;						
					modifiedMinorValences = HighCurvatureModification(spps2BeMerged, 'MINOR');
					seedPointsValence_(modifiedMinorValences,3) = 1;							
				end
			end
		end
	end

	numMinorPSLs = length(minorPSLpool_);
	for ii=1:numMinorPSLs
		minorPSL = minorPSLpool_(ii);
		if minorPSL.length>0	
			sppsEmptyMinorValence = find(0==seedPointsValence_(:,3));
            if ~isempty(sppsEmptyMinorValence)
				[potentialDisListMinor, potentialPosListMinor] = GetDisListOfPointList2Curve(...	
					seedPoints_(sppsEmptyMinorValence,:), [minorPSL.eleIndexList minorPSL.phyCoordList], 'MINOR');
				potentialSolidSppsMinor = find(potentialDisListMinor<=multiMergingThresholdsCtrl_(3)*relaxedFactor_);
				if ~isempty(potentialSolidSppsMinor)
					spps2BeMerged = sppsEmptyMinorValence(potentialSolidSppsMinor);
					seedPoints_(spps2BeMerged,:) = potentialPosListMinor(potentialSolidSppsMinor,:);
					seedPointsValence_(spps2BeMerged,3) = 1;
					modifiedMajorValences = HighCurvatureModification(spps2BeMerged, 'MAJOR');
					seedPointsValence_(modifiedMajorValences,1) = 1;	
					modifiedMediumValences = HighCurvatureModification(spps2BeMerged, 'MEDIUM');
					seedPointsValence_(modifiedMediumValences,2) = 1;						
				end
			end
		end
	end		
end

function modifiedValences = HighCurvatureModification(spps2BeMerged, psDir)
	global majorCoordList_;
    global mediumCoordList_;
	global minorCoordList_;		
	global seedPoints_;
	global seedPointsValence_;
	global mergingThreshold_;
	global relaxedFactor_;

	coordList = [];
	switch psDir
		case 'MAJOR'
			if isempty(majorCoordList_), modifiedValences = []; return; end
			coordList = majorCoordList_;
            spps2BeMerged = spps2BeMerged(0==seedPointsValence_(spps2BeMerged,1));
		case 'MEDIUM'
			if isempty(mediumCoordList_), modifiedValences = []; return; end
			coordList = mediumCoordList_;
            spps2BeMerged = spps2BeMerged(0==seedPointsValence_(spps2BeMerged,2));			
		case 'MINOR'
			if isempty(minorCoordList_), modifiedValences = []; return; end
			coordList = minorCoordList_;
            spps2BeMerged = spps2BeMerged(0==seedPointsValence_(spps2BeMerged,3));
	end
	pointList = seedPoints_(spps2BeMerged,end-2:end);
	disT = (coordList(:,1) - pointList(:,1)').^2;
	disT = disT + (coordList(:,2) - pointList(:,2)').^2;
	disT = disT + (coordList(:,3) - pointList(:,3)').^2;
	disT = sqrt(disT);		
	minVal = min(disT, [], 1);
	minVal = minVal/mergingThreshold_;
	modifiedValences = find(minVal<relaxedFactor_);	
	modifiedValences = spps2BeMerged(modifiedValences);
end

function [potentialDisList, potentialPosList] = GetDisListOfPointList2Curve(pointList, curveLine, psDir)
	global mergingThresholdMap_;
	disT = (curveLine(:,end-2) - pointList(:,end-2)').^2;
	disT = disT + (curveLine(:,end-1) - pointList(:,end-1)').^2;
	disT = disT + (curveLine(:,end) - pointList(:,end)').^2;
	disT = sqrt(disT);	
	[minVal, minValPos] = min(disT,[],1);

	potentialPosList = curveLine(minValPos,:);
	switch psDir
		case 'MAJOR', idx = 1;
		case 'MEDIUM', idx = 2;
		case 'MINOR', idx = 3;
	end
	
	potentialDisList = minVal';
	% potentialDisList = potentialDisList/mergingThreshold_;
	potentialDisList = potentialDisList ./ mergingThresholdMap_(potentialPosList(:,1),idx);
end

function [eleIndex, cartesianStress, principalStress, opt] = PreparingForTracing(startPoint)
	global nodeCoords_; 
	global eNodMat_;
	global cartesianStressField_;
	global eleCentroidList_;

	eleIndex = 0;
	cartesianStress = 0;
	principalStress = 0;	
	switch numel(startPoint)
		case 3
			disList = vecnorm(startPoint-eleCentroidList_, 2, 2);
			[~, targetEleIndex0] = min(disList);	
			[eleIndex, opt] = PositioningOnUnstructuredMesh(targetEleIndex0, startPoint);
			if ~opt, return; end			
		case 4
			eleIndex = startPoint(1);
			startPoint = startPoint(2:4);
            opt = 1;			
		otherwise
			error('Wrong Input For the Seed!')
	end
	NIdx = eNodMat_(eleIndex,:)';
	eleNodeCoords = nodeCoords_(NIdx,:);
	eleCartesianStress = cartesianStressField_(NIdx,:);			
	cartesianStress = ElementInterpolationInverseDistanceWeighting(eleNodeCoords, eleCartesianStress, startPoint);	
	principalStress = ComputePrincipalStress(cartesianStress);		
end

function iPSL = Have1morePSL(startPoint, tracingType)
	global integrationStepLimit_;

	iPSL = PrincipalStressLineStruct();
	switch tracingType
		case 'MAJOR', psDir = [10 11 12];
		case 'MEDIUM', psDir = [6 7 8];
		case 'MINOR', psDir = [2 3 4];
	end

	%%1. prepare for tracing			
	[eleIndex, cartesianStress, principalStress, opt] = PreparingForTracing(startPoint);
	if 0==opt, return; end
	
	%%2. tracing PSL
	startPoint = startPoint(end-2:end);
	PSLphyCoordList = startPoint;
	PSLcartesianStressList = cartesianStress;
	PSLeleIndexList = eleIndex;
	PSLprincipalStressList = principalStress;
	
	%%2.1 along first direction (v1)		
	[phyCoordList, cartesianStressList, eleIndexList, principalStressList] = ...
		TracingPSL_RK2(startPoint, principalStress(1,psDir), eleIndex, psDir, integrationStepLimit_);		
	PSLphyCoordList = [PSLphyCoordList; phyCoordList];
	PSLcartesianStressList = [PSLcartesianStressList; cartesianStressList];
	PSLeleIndexList = [PSLeleIndexList; eleIndexList];
	PSLprincipalStressList = [PSLprincipalStressList; principalStressList];
	
	%%2.2 along second direction (-v1)	
	[phyCoordList, cartesianStressList, eleIndexList, principalStressList] = ...
		TracingPSL_RK2(startPoint, -principalStress(1,psDir), eleIndex, psDir, integrationStepLimit_);		
	if size(phyCoordList,1) > 1
		phyCoordList = flip(phyCoordList);
		cartesianStressList = flip(cartesianStressList);
		eleIndexList = flip(eleIndexList);
		principalStressList = flip(principalStressList);
	end						
	PSLphyCoordList = [phyCoordList; PSLphyCoordList];
	PSLcartesianStressList = [cartesianStressList; PSLcartesianStressList];
	PSLeleIndexList = [eleIndexList; PSLeleIndexList];
	PSLprincipalStressList = [principalStressList; PSLprincipalStressList];
	
	%%2.3 finish Tracing the current major PSL	
	iPSL.midPointPosition = size(phyCoordList,1)+1;
	iPSL.length = size(PSLphyCoordList,1);
	iPSL.eleIndexList = PSLeleIndexList;
	iPSL.phyCoordList = PSLphyCoordList;
	iPSL.cartesianStressList = PSLcartesianStressList;	
	iPSL.principalStressList = PSLprincipalStressList;		
end

function val = PrincipalStressLineStruct()
	val = struct(...
		'ith',						0, 	...
		'length',					0,	...
		'midPointPosition',			0,	...		
		'phyCoordList',				[], ...
		'eleIndexList',				[], ...
		'cartesianStressList',		[],	...
		'principalStressList',		[] ...
	);	
end

function [phyCoordList, cartesianStressList, eleIndexList, principalStressList] = ...
			TracingPSL_RK2(startPoint, iniDir, elementIndex, typePSL, limiSteps)			
	global eNodMat_;
	global nodeCoords_;
	global cartesianStressField_;
	global tracingStepWidth_;
	phyCoordList = zeros(limiSteps,3);
	cartesianStressList = zeros(limiSteps,6);
	eleIndexList = zeros(limiSteps,1);
	principalStressList = zeros(limiSteps,12);	
	
	index = 0;
	k1 = iniDir;
	%%re-scale stepsize if necessary
	stepsize = tracingStepWidth_(elementIndex);
	testPot = startPoint + k1*tracingStepWidth_(elementIndex);
	[~, testPot1, bool1] = SearchNextIntegratingPointOnUnstructuredMesh(elementIndex, testPot, startPoint, 1);
	if bool1
		stepsize = norm(testPot1-startPoint)/norm(testPot-startPoint) * stepsize;
		%%initialize initial k1 and k2
		midPot = startPoint + k1*stepsize/2;
		[elementIndex2, ~, bool2] = SearchNextIntegratingPointOnUnstructuredMesh(elementIndex, midPot, startPoint, 0);
		if bool2 %%just in case
			NIdx = eNodMat_(elementIndex2,:)';
			vtxStress = cartesianStressField_(NIdx, :);
			vtxCoords = nodeCoords_(NIdx,:);
			cartesianStressOnGivenPoint = ElementInterpolationInverseDistanceWeighting(vtxCoords, vtxStress, midPot);
			principalStress = ComputePrincipalStress(cartesianStressOnGivenPoint);
			[k2, ~] = BidirectionalFeatureProcessing(k1, principalStress(typePSL));
			nextPoint = startPoint + stepsize*k2;
			[elementIndex, ~, bool3] = SearchNextIntegratingPointOnUnstructuredMesh(elementIndex, nextPoint, startPoint, 0);
			while bool3
				index = index + 1; if index > limiSteps, index = index-1; break; end
				NIdx = eNodMat_(elementIndex,:)';
				vtxStress = cartesianStressField_(NIdx, :);
				vtxCoords = nodeCoords_(NIdx,:); 
				cartesianStressOnGivenPoint = ElementInterpolationInverseDistanceWeighting(vtxCoords, vtxStress, nextPoint); 
				principalStress = ComputePrincipalStress(cartesianStressOnGivenPoint);	
				%%k1
				[k1, terminationCond] = BidirectionalFeatureProcessing(iniDir, principalStress(typePSL));	
				if ~terminationCond, index = index-1; break; end
				%%k2
				%%re-scale stepsize if necessary
				stepsize = tracingStepWidth_(elementIndex);
				testPot = nextPoint + k1*stepsize;
				[~, testPot1, bool1] = SearchNextIntegratingPointOnUnstructuredMesh(elementIndex, testPot, nextPoint, 1);			
				if ~bool1, index = index-1; break; end
				stepsize = norm(testPot1-nextPoint)/norm(testPot-nextPoint) * stepsize;
				midPot = nextPoint + k1*stepsize/2;
				[elementIndex2, ~, bool1] = SearchNextIntegratingPointOnUnstructuredMesh(elementIndex, midPot, nextPoint, 0);					
				if ~bool1, index = index-1; break; end	
				NIdx2 = eNodMat_(elementIndex2,:)';	
				vtxStress2 = cartesianStressField_(NIdx2,:);
				vtxCoords2 = nodeCoords_(NIdx2,:);
				cartesianStressOnGivenPoint2 = ElementInterpolationInverseDistanceWeighting(vtxCoords2, vtxStress2, midPot);
				principalStress2 = ComputePrincipalStress(cartesianStressOnGivenPoint2);
				[k2, ~] = BidirectionalFeatureProcessing(k1, principalStress2(typePSL));					
				%%store	
				iniDir = k1;
				phyCoordList(index,:) = nextPoint;
				cartesianStressList(index,:) = cartesianStressOnGivenPoint;
				eleIndexList(index,:) = elementIndex;
				principalStressList(index,:) = principalStress;
				%%next point
				nextPoint0 = nextPoint + stepsize*k2;
				[elementIndex, ~, bool3] = SearchNextIntegratingPointOnUnstructuredMesh(elementIndex, nextPoint0, nextPoint, 0);			
				nextPoint = nextPoint0;				
			end
		end
	end
	phyCoordList = phyCoordList(1:index,:);
	cartesianStressList = cartesianStressList(1:index,:);
	eleIndexList = eleIndexList(1:index,:);
	principalStressList = principalStressList(1:index,:);	
end

function [targetDirection, terminationCond] = BidirectionalFeatureProcessing(originalVec, Vec)
	global permittedMaxAdjacentTangentAngleDeviation_;
	terminationCond = 1;
	angle1 = acos(originalVec*Vec');
	angle2 = acos(-originalVec*Vec');	
	if angle1 < angle2
		targetDirection = Vec;
		if angle1 > pi/permittedMaxAdjacentTangentAngleDeviation_, terminationCond = 0; end
	else
		targetDirection = -Vec;
		if angle2 > pi/permittedMaxAdjacentTangentAngleDeviation_, terminationCond = 0; end
	end	
end

function [eleIndex, opt] = PositioningOnUnstructuredMesh(targetEleIndex0, startPoint)
	global eNodMat_; 
	global nodStruct_;
	global eleCentroidList_;
	opt = IsThisPointWithinThatElement(targetEleIndex0, startPoint);
	if opt
		eleIndex = targetEleIndex0;		
	else %% Search the Adjacent Elements
		tarNodes = eNodMat_(targetEleIndex0,:);
		allPotentialAdjacentElements = unique([nodStruct_(tarNodes(:)).adjacentEles]);
		potentialAdjacentElements = setdiff(allPotentialAdjacentElements, targetEleIndex0);
		for ii=1:length(potentialAdjacentElements)
			iEle = potentialAdjacentElements(ii);
			opt = IsThisPointWithinThatElement(iEle, startPoint);
			if opt, eleIndex = iEle; break; end
		end
	end
	if 0==opt		
		disList = vecnorm(startPoint-eleCentroidList_(allPotentialAdjacentElements,:), 2, 2);
		[~, nearOptimalEle] = min(disList);
		eleIndex = allPotentialAdjacentElements(nearOptimalEle);
	end
end

function opt = IsThisPointWithinThatElement(tarEleIndex, iCoord)
	global eleStruct_;
	global nodeCoords_; 
	global eNodMat_;
	global numFacesPerEle_;
	opt = 1; 
	iElefaceCentres = eleStruct_(tarEleIndex).faceCentres;
	iNodeCords = nodeCoords_(eNodMat_(tarEleIndex,:),:);
	%%compute direction vectors from iCoord to face centers as reference vectors
	refVec = iElefaceCentres - iCoord; %% dir vecs from volume center to face centers
	refVec2Vertices = iNodeCords - iCoord; 
	refVecNorm = vecnorm(refVec,2,2);
	refVec2VerticesNorm = vecnorm(refVec2Vertices,2,2);
	if 6==numFacesPerEle_
		thresholdAng = 91.0;
	else
		thresholdAng = 90.0;
	end
	if isempty(find(0==[refVecNorm; refVec2VerticesNorm])) %% iCoord does NOT coincides with a vertex or face center
		refVec = refVec ./ refVecNorm; 
		normVecs = eleStruct_(tarEleIndex).faceNormals;
		
		%%compute angle deviation
		angleDevs = zeros(numFacesPerEle_,1);
		for ii=1:numFacesPerEle_
			angleDevs(ii) = acos(refVec(ii,:)*normVecs(ii,:)')/pi*180;
		end
		%% iCoord is out of tarEleIndex, using the relaxed 91 instead of 90 for numerical instability
		maxAngle = max(angleDevs);
		if maxAngle > thresholdAng, opt = 0; end	
	end 
end

function val = ElementInterpolationInverseDistanceWeighting(coords, vtxEntity, ips)
	%% Inverse Distance Weighting
	%% coords --> element vertex coordinates, Matrix: [N-by-3] 
	%% vtxEntity --> entities on element vertics, Matrix: [N-by-M], e.g., M = 6 for 3D stress tensor
	%% ips --> to-be interpolated coordinate, Vector: [1-by-3]
	
	e = -2;
	D = vecnorm(ips-coords,2,2);
	[sortedD, sortedMapVec] = sort(D);
    if 0==sortedD(1)
        val = vtxEntity(sortedMapVec(1),:); return;
    end
	sortedVtxVals = vtxEntity(sortedMapVec,:);
	wV = sortedD.^e;
	V = sortedVtxVals.*wV;	
	val = sum(V) / sum(wV);
end

function principalStress = ComputePrincipalStress(cartesianStress)
	%% "cartesianStress" is in the order: Sigma_xx, Sigma_yy, Sigma_zz, Sigma_yz, Sigma_zx, Sigma_xy
	principalStress = zeros(1, 12);
	A = cartesianStress([1 6 5; 6 2 4; 5 4 3]);
	[eigenVec, eigenVal] = eig(A);
	principalStress([1 5 9]) = diag(eigenVal);
	principalStress([2 3 4 6 7 8 10 11 12]) = reshape(eigenVec,1,9);
end

function oPSLs = CompactStreamlines(iPSLs, truncatedThreshold)
	tarIndice = [];
	for ii=1:length(iPSLs)
		if iPSLs(ii).length > truncatedThreshold
			tarIndice(end+1,1) = ii;
		end
	end
	oPSLs = iPSLs(tarIndice);
	if isempty(oPSLs), oPSLs = []; end
end

function [nextElementIndex, p1, opt] = SearchNextIntegratingPointOnUnstructuredMesh(oldElementIndex, physicalCoordinates, sPoint, relocatingP1)
	global eleCentroidList_; 
	global eNodMat_; 
	global nodStruct_; 
	global eleState_;

	p1 = physicalCoordinates;
	nextElementIndex = oldElementIndex;	
	opt = IsThisPointWithinThatElement(oldElementIndex, p1);

	if opt
		return;
	else	
		tarNodes = eNodMat_(oldElementIndex,:); 
		potentialElements = unique([nodStruct_(tarNodes(:)).adjacentEles]);
		adjEleCtrs = eleCentroidList_(potentialElements,:);
		disList = vecnorm(p1-adjEleCtrs, 2, 2);
		[~, reSortMap] = sort(disList);
		potentialElements = potentialElements(reSortMap);
		for jj=1:length(potentialElements)
			iEle = potentialElements(jj);
			opt = IsThisPointWithinThatElement(iEle, p1);
			if opt, nextElementIndex = iEle; return; end
		end		
	end
	%%Scaling down the stepsize via Dichotomy
	if relocatingP1 && 0==opt	
		nn = 5;
		ii = 1;	
		while ii<=nn
			p1 = (sPoint+p1)/2;
			disList = vecnorm(p1-adjEleCtrs, 2, 2);
			[~, reSortMap] = sort(disList);
			potentialElements = potentialElements(reSortMap);			
			for jj=1:length(potentialElements)
				iEle = potentialElements(jj);
				opt = IsThisPointWithinThatElement(iEle, p1);
				if opt, nextElementIndex = iEle; return; end
			end
			ii = ii + 1;
		end
		if 0==eleState_(oldElementIndex)
			nextElementIndex = oldElementIndex; opt = 1;
		end
	end	
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Data Preparation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ImportStressFields(fileName)
	global boundingBox_;
	global meshType_;
	global numNodes_;
	global nodeCoords_;
	global numEles_;
	global eNodMat_;
	global numNodesPerEle_;
	global eleState_;
	global nodState_;
	global cartesianStressField_;
	global loadingCond_; 
	global fixingCond_;	
	global eleCentroidList_;
	global silhouetteStruct_;
	global eleSizeList_;
	
	global nodStruct_; 
	global eleStruct_; 
	global boundaryElements_; 
	global boundaryNodes_;
	global numNodesPerEle_;
	global numFacesPerEle_;
	global numNodesPerFace_;
	global eleFaceIndices_;	
	%%Read mesh and cartesian stress field
	fid = fopen(fileName, 'r');
	%%Mesh
	fgetl(fid); 
	domainType = fscanf(fid, '%s', 1);
	if ~strcmp(domainType, 'Solid'), warning('Un-supported Data!'); return; end
	meshType_ = fscanf(fid, '%s', 1);
	if ~(strcmp(meshType_, 'Hex') || strcmp(meshType_, 'Tet')), warning('Un-supported Mesh!'); return; end
	meshOrder = fscanf(fid, '%d', 1);
	if 1~=meshOrder, warning('Un-supported Mesh!'); return; end
	startReadingVertices = fscanf(fid, '%s', 1);
	if ~strcmp(startReadingVertices, 'Vertices:'), warning('Un-supported Data!'); return; end
	numNodes_ = fscanf(fid, '%d', 1);
	nodeCoords_ = fscanf(fid, '%e %e %e', [3, numNodes_])'; 
	startReadingElements = fscanf(fid, '%s', 1);
	if ~strcmp(startReadingElements, 'Elements:'), warning('Un-supported Data!'); return; end
	numEles_ = fscanf(fid, '%d', 1);
	switch meshType_
		case 'Hex'
			eNodMat_ = fscanf(fid, '%d %d %d %d %d %d %d %d', [8, numEles_])'; 
		case 'Tet'
			eNodMat_ = fscanf(fid, '%d %d %d %d', [4, numEles_])'; 
	end

	startReadingLoads = fscanf(fid, '%s %s', 2); 
	if ~strcmp(startReadingLoads, 'NodeForces:'), warning('Un-supported Data!'); return; end
	numLoadedNodes = fscanf(fid, '%d', 1);
	if numLoadedNodes>0, loadingCond_ = fscanf(fid, '%d %e %e %e', [4, numLoadedNodes])'; else, loadingCond_ = []; end
    
	startReadingFixations = fscanf(fid, '%s %s', 2);
    if ~strcmp(startReadingFixations, 'FixedNodes:'), warning('Un-supported Data!'); return; end
	numFixedNodes = fscanf(fid, '%d', 1);
	if numFixedNodes>0, fixingCond_ = fscanf(fid, '%d', [1, numFixedNodes])'; else, fixingCond_ = []; end
    
	startReadingStress = fscanf(fid, '%s %s', 2); 
	if ~strcmp(startReadingStress, 'CartesianStress:'), warning('Un-supported Data!'); return; end
	numValidNods = fscanf(fid, '%d', 1);
	cartesianStressField_ = fscanf(fid, '%e %e %e %e %e %e', [6, numValidNods])';		
	fclose(fid);
		
	%%Extract Boundary Element Info.
	switch meshType_
		case 'Hex'
			numNodesPerEle_ = 8;
			numFacesPerEle_ = 6;
			numNodesPerFace_ = 4;
			eleFaceIndices_ = [4 3 2 1; 5 6 7 8; 1 2 6 5; 8 7 3 4; 5 8 4 1; 2 3 7 6];
			pp = [1 2 3 4];	
		case 'Tet'
			numNodesPerEle_ = 4;
			numFacesPerEle_ = 4;
			numNodesPerFace_ = 3;
			eleFaceIndices_ = [3 2 1;  1 2 4;  2 3 4;  3 1 4];
			pp = [1 2 2 3];	
	end
	
	boundingBox_ = [min(nodeCoords_, [], 1); max(nodeCoords_, [], 1)];	
	[surfaceMeshElements_, surfaceMeshNodeCoords_, nodState_, boundaryNodes_] = ExtractBoundaryInfoFromSolidMesh();
	%%Extract Silhouette for Vis.
	silhouetteStruct_.vertices = surfaceMeshNodeCoords_;
	silhouetteStruct_.faces = surfaceMeshElements_;

	%%element centroids
	eleNodCoordListX = nodeCoords_(:,1); eleNodCoordListX = eleNodCoordListX(eNodMat_);
	eleNodCoordListY = nodeCoords_(:,2); eleNodCoordListY = eleNodCoordListY(eNodMat_);
	eleNodCoordListZ = nodeCoords_(:,3); eleNodCoordListZ = eleNodCoordListZ(eNodMat_);
	eleCentroidList_ = [sum(eleNodCoordListX,2) sum(eleNodCoordListY,2) sum(eleNodCoordListZ,2)]/numNodesPerEle_;		
	
	%% Build Element Tree for Unstructured Quad-Mesh	
	iNodStruct = struct('adjacentEles', []); 
	nodStruct_ = repmat(iNodStruct, numNodes_, 1);
	for ii=1:numEles_
		for jj=1:numNodesPerEle_
			nodStruct_(eNodMat_(ii,jj)).adjacentEles(1,end+1) = ii;
		end
	end		
	boundaryElements_ = unique([nodStruct_(boundaryNodes_).adjacentEles]);
	boundaryElements_ = boundaryElements_(:);
	eleState_ = zeros(numEles_,1);
	eleState_(boundaryElements_,1) = 1;		
	%% build element tree		
	iEleStruct = struct('faceCentres', [], 'faceNormals', []); %%pure-Hex
	eleStruct_ = repmat(iEleStruct, numEles_, 1);
	for ii=1:numEles_
		iNodes = eNodMat_(ii,:);
		iEleVertices = nodeCoords_(iNodes, :);
		iEleFacesX = iEleVertices(:,1); iEleFacesX = iEleFacesX(eleFaceIndices_);
		iEleFacesY = iEleVertices(:,2); iEleFacesY = iEleFacesY(eleFaceIndices_);
		iEleFacesZ = iEleVertices(:,3); iEleFacesZ = iEleFacesZ(eleFaceIndices_);
		ACs = [iEleFacesX(:,pp(1))-iEleFacesX(:,pp(3)) iEleFacesY(:,pp(1))-iEleFacesY(:,pp(3)) iEleFacesZ(:,pp(1))-iEleFacesZ(:,pp(3))];
		BDs = [iEleFacesX(:,pp(2))-iEleFacesX(:,pp(4)) iEleFacesY(:,pp(2))-iEleFacesY(:,pp(4)) iEleFacesZ(:,pp(2))-iEleFacesZ(:,pp(4))];
		iACxBD = cross(ACs,BDs); 
		aveNormal = iACxBD ./ vecnorm(iACxBD,2,2);			
		tmp = iEleStruct;			
		%% tmp.faceNormals = aveNormal;
		%% in case the node orderings on each element face are not constant
		tmp.faceCentres = [sum(iEleFacesX,2) sum(iEleFacesY,2) sum(iEleFacesZ,2)]/numNodesPerFace_;
		iEleCt = eleCentroidList_(ii,:);
		refVecs = iEleCt - tmp.faceCentres; refVecs = refVecs ./ vecnorm(refVecs,2,2);
		dirEval = acos(sum(refVecs .* aveNormal, 2));
		dirDes = ones(numFacesPerEle_,1); dirDes(dirEval<pi/2) = -1;
		faceNormals = dirDes .* aveNormal;
		tmp.faceNormals = faceNormals;
		eleStruct_(ii) = tmp;
	end
	
	%% Evaluate Element Sizes
	tmpSizeList = zeros(numFacesPerEle_, numEles_);
	for ii=1:numEles_
		tmpSizeList(:,ii) = vecnorm(eleCentroidList_(ii,:)-eleStruct_(ii).faceCentres,2,2);
	end
	eleSizeList_ = 2*min(tmpSizeList,[],1)';		
end

function [boundaryFaceNodMat, boundaryFaceNodeCoords, nodState, boundaryNodes] = ExtractBoundaryInfoFromSolidMesh()
	global numNodes_;
	global numEles_;
	global eNodMat_;
	global nodeCoords_;
	global numFacesPerEle_;
	global numNodesPerFace_;
	global eleFaceIndices_;

	eleFaces = eleFaceIndices_'; eleFaces = eleFaces(:)';
	patchIndices = eNodMat_(:,eleFaces)';
	patchIndices = reshape(patchIndices(:), numNodesPerFace_, numFacesPerEle_*numEles_)';	
	tmp = sort(patchIndices,2);
	[uniqueFaces, ia, ~] = unique(tmp, 'stable', 'rows');
	leftFaceIDs = (1:numFacesPerEle_*numEles_)'; leftFaceIDs = setdiff(leftFaceIDs, ia);
	leftFaces = tmp(leftFaceIDs,:);
	[~, surfFacesIDsInUniqueFaces] = setdiff(uniqueFaces, leftFaces, 'rows');
	boundaryFaceNodMat = patchIndices(ia(surfFacesIDsInUniqueFaces),:);
	boundaryNodes = int32(unique(boundaryFaceNodMat));
	nodState = zeros(numNodes_,1); nodState(boundaryNodes) = 1;	
	
	allNodes = zeros(numNodes_,1);
	allNodes(boundaryNodes) = (1:numel(boundaryNodes))';
	boundaryFaceNodMat = allNodes(boundaryFaceNodMat);
	boundaryFaceNodeCoords = nodeCoords_(boundaryNodes,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowPSLs(varargin)
    global boundingBox_;
    global silhouetteStruct_;
	global majorPSLpool_;
	global mediumPSLpool_;
	global minorPSLpool_;
	
	numTarMajorPSLs = numel(majorPSLpool_);
	numTarMediumPSLs = numel(mediumPSLpool_);
	numTarMinorPSLs = numel(minorPSLpool_);
	color4MajorPSLs = struct('arr', []); color4MajorPSLs = repmat(color4MajorPSLs, numTarMajorPSLs, 1);
	color4MediumPSLs = struct('arr', []); color4MediumPSLs = repmat(color4MediumPSLs, numTarMediumPSLs, 1);
	color4MinorPSLs = struct('arr', []); color4MinorPSLs = repmat(color4MinorPSLs, numTarMinorPSLs, 1);	
	for ii=1:numTarMajorPSLs
		color4MajorPSLs(ii).arr = ones(1, majorPSLpool_(ii).length);
	end
	for ii=1:numTarMediumPSLs
		color4MediumPSLs(ii).arr = ones(1, mediumPSLpool_(ii).length);
	end			
	for ii=1:numTarMinorPSLs
		color4MinorPSLs(ii).arr = ones(1, minorPSLpool_(ii).length);
	end	
	if 0==nargin, thicknessScaling = 150; else, thicknessScaling = varargin{1}; end
	lineWidthTube = min(boundingBox_(2,:)-boundingBox_(1,:))/thicknessScaling;
	[gridXmajor, gridYmajor, gridZmajor, gridCmajor, ~] = ExpandPSLs2Tubes(majorPSLpool_, color4MajorPSLs, lineWidthTube);
	[gridXmedium, gridYmedium, gridZmedium, gridCmedium, ~] = ExpandPSLs2Tubes(mediumPSLpool_, color4MediumPSLs, lineWidthTube);
	[gridXminor, gridYminor, gridZminor, gridCminor, ~] = ExpandPSLs2Tubes(minorPSLpool_, color4MinorPSLs, lineWidthTube);
	handleMajorPSL = []; handleMediumPSL = []; handleMinorPSL = [];
	%%Show silhouette
	hSilo = patch(gca, silhouetteStruct_); hold(gca, 'on');	
	if ~isempty(gridXmajor)
		hold(gca, 'on'); 
		handleMajorPSL = surf(gca, gridXmajor, gridYmajor, gridZmajor, gridCmajor);
	end
	if ~isempty(gridXmedium)
		hold(gca, 'on');
		handleMediumPSL = surf(gca, gridXmedium, gridYmedium, gridZmedium, gridCmedium);
	end
	if ~isempty(gridXminor)
		hold(gca, 'on');
		handleMinorPSL = surf(gca, gridXminor, gridYminor, gridZminor, gridCminor);
	end	
	set(handleMajorPSL, 'FaceColor', [252 141 98]/255, 'EdgeColor', 'None');
	set(handleMediumPSL, 'FaceColor', [0 176 80]/255, 'EdgeColor', 'None');
	set(handleMinorPSL, 'FaceColor', [91 155 213]/255, 'EdgeColor', 'None');
	set(hSilo, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.1, 'EdgeColor', 'None');
	view(gca, 3);
	camproj(gca, 'perspective');
	axis(gca, 'equal'); 
	axis(gca, 'tight');
	axis(gca, 'off');
	%%Lighting, Reflection
	lighting(gca, 'gouraud');
	material(gca, 'shiny');
	camlight(gca, 'headlight', 'infinite');    
end

function [gridX, gridY, gridZ, gridC, gridIndices] = ExpandPSLs2Tubes(PSLs, colorSrc, r)
	%%Syntax
	%% [gridX, gridY, gridZ, gridC] = ExpandPSLs2Tubes(PSLs, colorSrc, r)
	gridX = [];
	gridY = [];
	gridZ = [];
	gridC = [];
	gridIndices = [];
	if isempty(PSLs), return; end
	n = 8; 
	numLines = length(PSLs);
	gridXYZ = zeros(3,n+1,1);
	gridC = zeros(n+1,1);
	gridIndices = struct('arr', []);
	gridIndices = repmat(gridIndices, numLines, 1);	
	for ii=1:numLines		
		curve = PSLs(ii).phyCoordList';
		npoints = size(curve,2);
		%deltavecs: average for internal points. first strecth for endpoitns.		
		dv = curve(:,[2:end,end])-curve(:,[1,1:end-1]);		
		%make nvec not parallel to dv(:,1)
		nvec=zeros(3,1); 
		[~,idx]=min(abs(dv(:,1))); 
		nvec(idx)=1;
		%precalculate cos and sing factors:
		cfact=repmat(cos(linspace(0,2*pi,n+1)),[3,1]);
		sfact=repmat(sin(linspace(0,2*pi,n+1)),[3,1]);
		%Main loop: propagate the normal (nvec) along the tube
		xyz = zeros(3,n+1,npoints+2);
		for k=1:npoints
			convec=cross(nvec,dv(:,k));
			convec=convec./norm(convec);
			nvec=cross(dv(:,k),convec);
			nvec=nvec./norm(nvec);
			%update xyz:
			xyz(:,:,k+1)=repmat(curve(:,k),[1,n+1]) + cfact.*repmat(r*nvec,[1,n+1]) + sfact.*repmat(r*convec,[1,n+1]);
        end
		%finally, cap the ends:
		xyz(:,:,1)=repmat(curve(:,1),[1,n+1]);
		xyz(:,:,end)=repmat(curve(:,end),[1,n+1]);
		gridIndices(ii).arr = size(gridXYZ,3) : size(gridXYZ,3)+size(xyz,3)-1;
		gridXYZ(:,:,end+1:end+npoints+2) = xyz;	
		color = colorSrc(ii).arr;	
		c = [color(1) color color(end)];
		c = repmat(c, n+1, 1);
		gridC(:,end+1:end+npoints+2) = c;
	end		
	gridX = squeeze(gridXYZ(1,:,:)); 
	gridX(:,1) = [];
	gridY = squeeze(gridXYZ(2,:,:)); 
	gridY(:,1) = [];
	gridZ = squeeze(gridXYZ(3,:,:)); 
	gridZ(:,1) = [];
	gridC(:,1) = [];
end

function ShowProblemDescription()
	global nodeCoords_;
	global loadingCond_;
	global fixingCond_;
	global boundingBox_;
	global silhouetteStruct_;

	
	hd = patch(silhouetteStruct_); hold('on');
	set(hd, 'FaceColor', [65 174 118]/255, 'FaceAlpha', 0.9, 'EdgeColor', 'k');
    if ~isempty(loadingCond_)
		lB = 0.2; uB = 1.0;
		amps = vecnorm(loadingCond_(:,2:end),2,2);
		maxAmp = max(amps); minAmp = min(amps);
		if abs(minAmp-maxAmp)/(minAmp+maxAmp)<0.1
			scalingFac = 1;
		else
			if minAmp/maxAmp>lB/uB, lB = minAmp/maxAmp; end
			scalingFac = lB + (uB-lB)*(amps-minAmp)/(maxAmp-minAmp);
		end
		loadingDirVec = loadingCond_(:,2:end)./amps.*scalingFac;
		coordLoadedNodes = nodeCoords_(loadingCond_(:,1),:);
		amplitudesF = mean(boundingBox_(2,:)-boundingBox_(1,:))/5 * loadingDirVec;
		hold('on'); quiver3(coordLoadedNodes(:,1), coordLoadedNodes(:,2), coordLoadedNodes(:,3), amplitudesF(:,1), ...
			amplitudesF(:,2), amplitudesF(:,3), 0, 'Color', [255 127 0.0]/255, 'LineWidth', 2, 'MaxHeadSize', 1, 'MaxHeadSize', 1);
	end
    if ~isempty(fixingCond_)
		tarNodeCoord = nodeCoords_(fixingCond_(:,1),:);
		hold('on'); hd1 = plot3(tarNodeCoord(:,1), tarNodeCoord(:,2), tarNodeCoord(:,3), 'x', ...
			'color', [153 153 153]/255, 'LineWidth', 3, 'MarkerSize', 15);		
	end
	view(gca, 3);
	camproj(gca, 'perspective');
	axis(gca, 'equal'); 
	axis(gca, 'tight');
	axis(gca, 'off');
	
	%%Lighting, Reflection
	lighting(gca, 'gouraud');
	material(gca, 'dull');
	camlight(gca, 'headlight', 'infinite');
end