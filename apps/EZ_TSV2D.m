%% EZ_TSV2D.m
%% This is an Easy-to-Use version of 3D-TSV for 2D cases, 
%% all of the functionalities are included in this single script.
%% It supports to visualize the stress tensor field simulated on the 1st-order Quadrilateral and Triangular Mesh.
%% A demo-dataset is provided in '../data/demo_data_2D.stress'
% Author: Junpeng Wang (junpeng.wang@tum.de)
% Date: 2023-12-12
clear; clc;

global majorPSLpool_; 
global minorPSLpool_; 
global degeneratePoints_;

%%Extral Ctrls
global TopologyAnalysis_;
global excludeDegeneratePointsOnBoundary_;

%%1. Import Data
stressfileName = '../data/demoData_2D_Tri_Femur.TSV';
ImportStressFields(stressfileName);
figure; ShowProblemDescription();

%%2. Additional Settings
TopologyAnalysis_ = 1; %% 1 == 'ON', 0 == 'OFF'
excludeDegeneratePointsOnBoundary_ = 0; 

%%3. PSLs generation
PSLsDensityCtrl = 20; %%The larger this value is, the denser the PSLs are.
TSV2D(PSLsDensityCtrl);
figure; ShowPSLs();


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TSV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TSV2D(PSLsDensityCtrl, varargin)
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
    global minorPSLpool_; 
	global majorCoordList_; 
    global minorCoordList_;
	global TopologyAnalysis_;
	global degeneratePoints_;
	global nodeWiseStressField_;
	
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
	seedPointsValence_ = zeros(size(seedPointsHistory_));
	numSeedPoints = size(seedPointsHistory_,1);
	seedPoints_ = [seedAssociatedEles_ seedPointsHistory_];
	
	InitializeMergingThresholdMap();
	
	majorPSLpool_ = PrincipalStressLineStruct();
	minorPSLpool_ = PrincipalStressLineStruct();
	degeneratePoints_ = [];
	if TopologyAnalysis_ && nodeWiseStressField_, TensorFieldTopologyAnalysis(); end
	PreprocessSeedPoints();
	
	majorCoordList_ = [];	
	minorCoordList_ = [];
	
	startCoord = boundingBox_(1,:) + sum(boundingBox_, 1)/2;
	
	%% Seeding
	its = 0;
	looper = sum(sum(seedPointsValence_));
	while looper<2*numSeedPoints
		its = its + 1;
		valenceMetric = sum(seedPointsValence_,2);
		unFinishedSpps = find(valenceMetric<2);
		spp = unFinishedSpps(1);		
		if 0==looper
			[~, spp] = min(vecnorm(startCoord-seedPoints_(:,end-1:end),2,2));
		else
			tmp	= seedPointsValence_(unFinishedSpps,:);
			tmp = find(1==sum(tmp,2));
			if ~isempty(tmp), unFinishedSpps = unFinishedSpps(tmp); end					
			[~, tarPos] = min(vecnorm(startCoord-seedPoints_(unFinishedSpps,end-1:end),2,2));					
			spp = unFinishedSpps(tarPos);		
		end		
		valences = seedPointsValence_(spp,:);						
		seed = seedPoints_(spp,:);
		if 0==valences(1)
			seedPointsValence_(spp,1) = 1;
			majorPSL = Have1morePSL2D(seed, 'MAJOR', []);		
			if 0==majorPSL.length
				looper = sum(sum(seedPointsValence_)); 
				disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
					' Total.: ' sprintf('%6i',2*numSeedPoints)]); continue; 
			end			
			majorPSLpool_(end+1,1) = majorPSL;				
			majorCoordList_(end+1:end+majorPSL.length,:) = majorPSL.phyCoordList;
			sppsEmptyMajorValence = find(0==seedPointsValence_(:,1));
			if ~isempty(sppsEmptyMajorValence)
				[potentialDisListMajor, potentialPosListMajor] = GetDisListOfPointList2Curve2D(seedPoints_(...
						sppsEmptyMajorValence,:), [majorPSL.eleIndexList majorPSL.phyCoordList], 'MAJOR');					
				potentialSolidSppsMajor = find(potentialDisListMajor<relaxedFactor_);
				if ~isempty(potentialSolidSppsMajor)
					spps2BeMerged = sppsEmptyMajorValence(potentialSolidSppsMajor);
					seedPoints_(spps2BeMerged,:) = potentialPosListMajor(potentialSolidSppsMajor,:);								
					seedPointsValence_(spps2BeMerged,1) = 1;				
					modifiedMinorValences = HighCurvatureModification2D(spps2BeMerged, 'MINOR');				
					seedPointsValence_(modifiedMinorValences,2) = 1;	
				end
			end				
		end
			
		if 0==valences(2)
			seedPointsValence_(spp,2) = 1;			
			minorPSL = Have1morePSL2D(seed, 'MINOR', []);	
			if 0==minorPSL.length
				looper = sum(sum(seedPointsValence_)); 
				disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
					' Total.: ' sprintf('%6i',2*numSeedPoints)]); continue; 
			end		
			minorPSLpool_(end+1,1) = minorPSL;
			minorCoordList_(end+1:end+minorPSL.length,:) = minorPSL.phyCoordList;				
			sppsEmptyMinorValence = find(0==seedPointsValence_(:,2));
			if ~isempty(sppsEmptyMinorValence)   
				[potentialDisListMinor, potentialPosListMinor] = GetDisListOfPointList2Curve2D(seedPoints_(...
						sppsEmptyMinorValence,:), [minorPSL.eleIndexList minorPSL.phyCoordList], 'MINOR');					
				potentialSolidSppsMinor = find(potentialDisListMinor<relaxedFactor_);
				if ~isempty(potentialSolidSppsMinor)
					spps2BeMerged = sppsEmptyMinorValence(potentialSolidSppsMinor);
					seedPoints_(spps2BeMerged,:) = potentialPosListMinor(potentialSolidSppsMinor,:);
					seedPointsValence_(spps2BeMerged,2) = 1;				
					modifiedMajorValences = HighCurvatureModification2D(spps2BeMerged, 'MAJOR');					
					seedPointsValence_(modifiedMajorValences,1) = 1;				
				end
			end					
		end
		looper = sum(sum(seedPointsValence_));
		disp([' Iteration.: ' sprintf('%4i',its) ' Progress.: ' sprintf('%6i',looper) ...
			' Total.: ' sprintf('%6i',2*numSeedPoints)]);			
	end

	majorPSLpool_ = CompactStreamlines(majorPSLpool_, 5);
	minorPSLpool_ = CompactStreamlines(minorPSLpool_, 5);	
end

function InitializeMergingThresholdMap()
	global numEles_;
	global mergingThreshold_;
	global mergingThresholdMap_;
	mergingThresholdMap_ = repmat(mergingThreshold_, numEles_, 2);
end

function PreprocessSeedPoints()
	global seedPoints_;
	global seedPointsValence_;
	global majorPSLpool_; 
    global minorPSLpool_; 
	global relaxedFactor_;
	numMajorPSLs = length(majorPSLpool_);
	for ii=1:numMajorPSLs
		majorPSL = majorPSLpool_(ii);
		if majorPSL.length>0					
			sppsEmptyMajorValence = find(0==seedPointsValence_(:,1));
            if ~isempty(sppsEmptyMajorValence)
				[potentialDisListMajor, potentialPosListMajor] = GetDisListOfPointList2Curve2D(...	
					seedPoints_(sppsEmptyMajorValence,:), [majorPSL.eleIndexList majorPSL.phyCoordList], 'MAJOR');
				potentialSolidSppsMajor = find(potentialDisListMajor<=relaxedFactor_);
				if ~isempty(potentialSolidSppsMajor)
					spps2BeMerged = sppsEmptyMajorValence(potentialSolidSppsMajor);							
					seedPoints_(spps2BeMerged,:) = potentialPosListMajor(potentialSolidSppsMajor,:);
					seedPointsValence_(spps2BeMerged,1) = 1;											
					modifiedMinorValences = HighCurvatureModification2D(spps2BeMerged, 'MINOR');
					seedPointsValence_(modifiedMinorValences,2) = 1;							
				end
			end
		end
	end
	numMinorPSLs = length(minorPSLpool_);
	for ii=1:numMinorPSLs
		minorPSL = minorPSLpool_(ii);
		if minorPSL.length>0	
			sppsEmptyMinorValence = find(0==seedPointsValence_(:,2));
            if ~isempty(sppsEmptyMinorValence)
				[potentialDisListMinor, potentialPosListMinor] = GetDisListOfPointList2Curve2D(...	
					seedPoints_(sppsEmptyMinorValence,:), [minorPSL.eleIndexList minorPSL.phyCoordList], 'MINOR');
				potentialSolidSppsMinor = find(potentialDisListMinor<=relaxedFactor_);
				if ~isempty(potentialSolidSppsMinor)
					spps2BeMerged = sppsEmptyMinorValence(potentialSolidSppsMinor);
					seedPoints_(spps2BeMerged,:) = potentialPosListMinor(potentialSolidSppsMinor,:);
					seedPointsValence_(spps2BeMerged,2) = 1;
					modifiedMajorValences = HighCurvatureModification2D(spps2BeMerged, 'MAJOR');
					seedPointsValence_(modifiedMajorValences,1) = 1;						
				end
			end
		end
	end
end

function modifiedValences = HighCurvatureModification2D(spps2BeMerged, psDir)
	global majorCoordList_; 
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
		case 'MINOR'
			if isempty(minorCoordList_), modifiedValences = []; return; end
			coordList = minorCoordList_;
            spps2BeMerged = spps2BeMerged(0==seedPointsValence_(spps2BeMerged,2));
	end
	pointList = seedPoints_(spps2BeMerged,end-1:end);
	disT = (coordList(:,1) - pointList(:,1)').^2;
	disT = disT + (coordList(:,2) - pointList(:,2)').^2;
	disT = sqrt(disT);		
	minVal = min(disT, [], 1);
	minVal = minVal/mergingThreshold_;
	modifiedValences = find(minVal<relaxedFactor_);	
	modifiedValences = spps2BeMerged(modifiedValences);
end

function [potentialDisList, potentialPosList] = GetDisListOfPointList2Curve2D(pointList, curveLine, psDir)
	global mergingThresholdMap_;	
	disT = (curveLine(:,end-1) - pointList(:,end-1)').^2;
	disT = disT + (curveLine(:,end) - pointList(:,end)').^2;
	disT = sqrt(disT);	
	[minVal, minValPos] = min(disT,[],1);
	
	potentialPosList = curveLine(minValPos,:);
	switch psDir
		case 'MAJOR', idx = 1;	
		case 'MINOR', idx = 2;
	end
	
	potentialDisList = minVal';
	% potentialDisList = potentialDisList/mergingThreshold_;
	potentialDisList = potentialDisList ./ mergingThresholdMap_(potentialPosList(:,1),idx);
end

function iPSL = Have1morePSL2D(startPoint, tracingType, iniDir)
	global integrationStepLimit_;

	iPSL = PrincipalStressLineStruct();
	switch tracingType
		case 'MAJOR', psDir = [5 6];	
		case 'MINOR', psDir = [2 3];
	end

	%%1. prepare for tracing			
	[eleIndex, cartesianStress, principalStress, opt] = PreparingForTracing2D(startPoint);
	if 0==opt, return; end
	
	if ~isempty(iniDir)
		principalStress(1,psDir) = iniDir;
	end
	
	%%2. tracing PSL
	startPoint = startPoint(end-1:end);
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
		'importanceMetric',			[],	...
		'principalStressList',		[] ...
	);	
end

function [phyCoordList, cartesianStressList, eleIndexList, principalStressList] = ...
			TracingPSL_RK2(startPoint, iniDir, elementIndex, typePSL, limiSteps)
	global eNodMat_;
	global nodeCoords_;
	global cartesianStressField_;
	global tracingStepWidth_;
	global nodeWiseStressField_;
	
	phyCoordList = zeros(limiSteps,2);
	cartesianStressList = zeros(limiSteps,3);
	eleIndexList = zeros(limiSteps,1);
	principalStressList = zeros(limiSteps,6);	
	
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
			if nodeWiseStressField_
				NIdx = eNodMat_(elementIndex2,:)';
				vtxStress = cartesianStressField_(NIdx, :);
				vtxCoords = nodeCoords_(NIdx,:);
				cartesianStressOnGivenPoint = ElementInterpolationInverseDistanceWeighting(vtxCoords, vtxStress, midPot);
			else
				cartesianStressOnGivenPoint = cartesianStressField_(elementIndex2, :); %% Constant in each element
			end
			principalStress = ComputePrincipalStress2D(cartesianStressOnGivenPoint);
			[k2, ~] = BidirectionalFeatureProcessing(k1, principalStress(typePSL));
			nextPoint = startPoint + stepsize*k2;
			[elementIndex, ~, bool3] = SearchNextIntegratingPointOnUnstructuredMesh(elementIndex, nextPoint, startPoint, 0);
			while bool3
				index = index + 1; if index > limiSteps, index = index-1; break; end
				if nodeWiseStressField_
					NIdx = eNodMat_(elementIndex,:)';
					vtxStress = cartesianStressField_(NIdx, :);
					vtxCoords = nodeCoords_(NIdx,:); 
					cartesianStressOnGivenPoint = ElementInterpolationInverseDistanceWeighting(vtxCoords, vtxStress, nextPoint); 				
				else
					cartesianStressOnGivenPoint = cartesianStressField_(elementIndex, :);
				end
				principalStress = ComputePrincipalStress2D(cartesianStressOnGivenPoint);
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
				if nodeWiseStressField_
					NIdx2 = eNodMat_(elementIndex2,:)';	
					vtxStress2 = cartesianStressField_(NIdx2,:);
					vtxCoords2 = nodeCoords_(NIdx2,:);
					cartesianStressOnGivenPoint2 = ElementInterpolationInverseDistanceWeighting(vtxCoords2, vtxStress2, midPot);				
				else
					cartesianStressOnGivenPoint2 = cartesianStressField_(elementIndex2,:);
				end
				principalStress2 = ComputePrincipalStress2D(cartesianStressOnGivenPoint2);
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

function [targetDirection, terminationCond] = BidirectionalFeatureProcessing(originalVec, Vecs)
	global permittedMaxAdjacentTangentAngleDeviation_;
	potentialVecs = [Vecs; -Vecs];
	angList = acos(potentialVecs*originalVec');
	[minAng, minAngPos] = min(angList);
	targetDirection = potentialVecs(minAngPos,:);
	if minAng < pi/permittedMaxAdjacentTangentAngleDeviation_
		terminationCond = 1;
	else
		terminationCond = 0;
	end
end

function [eleIndex, cartesianStress, principalStress, opt] = PreparingForTracing2D(startPoint)
	global nodeCoords_; 
	global eNodMat_; 
	global cartesianStressField_;
	global eleCentroidList_;
	global nodeWiseStressField_;
	
	eleIndex = 0;
	cartesianStress = 0;
	principalStress = 0;
	switch numel(startPoint)
		case 2
			disList = vecnorm(startPoint-eleCentroidList_, 2, 2);
			[~, targetEleIndex0] = min(disList);	
			[eleIndex, opt] = PositioningOnUnstructuredMesh(targetEleIndex0, startPoint);
			if ~opt, return; end			
		case 3
			eleIndex = startPoint(1);
			startPoint = startPoint(2:3);
            opt = 1;			
		otherwise
			error('Wrong Input For the Seed!')
	end
	if nodeWiseStressField_
		NIdx = eNodMat_(eleIndex,:)';
		eleNodeCoords = nodeCoords_(NIdx,:);
		eleCartesianStress = cartesianStressField_(NIdx,:);			
		cartesianStress = ElementInterpolationInverseDistanceWeighting(eleNodeCoords, eleCartesianStress, startPoint);	
	else
		cartesianStress = cartesianStressField_(eleIndex,:);
	end
	principalStress = ComputePrincipalStress2D(cartesianStress);
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
	global numEdgesPerEle_;
	
	opt = 1; 
	iEleEdgeCentres = eleStruct_(tarEleIndex).edgeCentres;
	iNodeCords = nodeCoords_(eNodMat_(tarEleIndex,:),:);
	%%compute direction vectors from iCoord to edge centers as reference vectors
	refVec = iEleEdgeCentres - iCoord; %% dir vecs from volume center to edge centers
	refVec2Vertices = iNodeCords - iCoord; 
	refVecNorm = vecnorm(refVec,2,2);
	refVec2VerticesNorm = vecnorm(refVec2Vertices,2,2);
	if isempty(find(0==[refVecNorm; refVec2VerticesNorm])) %% iCoord does NOT coincides with a vertex or face center
		refVec = refVec ./ refVecNorm; 
		normVecs = eleStruct_(tarEleIndex).edgeNormals;
		
		%%compute angle deviation
		angleDevs = zeros(numEdgesPerEle_,1);
		for ii=1:numEdgesPerEle_
			angleDevs(ii) = acos(refVec(ii,:)*normVecs(ii,:)')/pi*180;
		end
		%% iCoord is out of tarEleIndex, using the relaxed 91 instead of 90 for numerical instability
		maxAngle = max(angleDevs);
		if maxAngle > 90.0, opt = 0; end	
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

function principalStress = ComputePrincipalStress2D(cartesianStress)
	iStressTensor = cartesianStress(1,:);
	iStressTensor = iStressTensor([1 3; 3 2]);
	[eigenVec, eigenVal] = eig(iStressTensor);
	principalStress = [eigenVal(1,1); eigenVec(:,1); eigenVal(2,2); eigenVec(:,2)]';	
end

function val = ComputeVonMisesStress2D(carStress)
	%% "carStress" is in the order: Sigma_xx, Sigma_yy, Sigma_xy
	val = sqrt(carStress(:,1).^2 + carStress(:,2).^2 - carStress(:,1).*carStress(:,2) + 3*carStress(:,3).^2 );
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
	global numEdgesPerEle_;
	global eleState_;
	global nodState_;
	global boundaryEdgeNodMat_;
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
	global nodeWiseStressField_;
	%%Read mesh and cartesian stress field
	fid = fopen(fileName, 'r');
	%%Mesh
	fgetl(fid); 
	tmp = fscanf(fid, '%s %s %s', 3);
	type = fscanf(fid, '%s', 1);
	switch type
		case 'NODE', nodeWiseStressField_ = 1;
		case 'ELEMENT', nodeWiseStressField_ = 0; %%Element-wise Stress Data
		otherwise, error('Un-supported Stress Data!');
	end		
	domainType = fscanf(fid, '%s', 1);
	if ~strcmp(domainType, 'Plane'), error('Un-supported Data!'); end
	meshType_ = fscanf(fid, '%s', 1);
	if ~(strcmp(meshType_, 'Quad') || strcmp(meshType_, 'Tri')), error('Un-supported Mesh!'); end
	meshOrder = fscanf(fid, '%d', 1);
	if 1~=meshOrder, warning('Un-supported Mesh!'); return; end
	startReadingVertices = fscanf(fid, '%s', 1);
	if ~strcmp(startReadingVertices, 'Vertices:'), error('Un-supported Data!'); end
	numNodes_ = fscanf(fid, '%d', 1);
	nodeCoords_ = fscanf(fid, '%e %e', [2, numNodes_])'; 
	startReadingElements = fscanf(fid, '%s', 1);
	if ~strcmp(startReadingElements, 'Elements:'), error('Un-supported Data!'); end
	numEles_ = fscanf(fid, '%d', 1);
	switch meshType_
		case 'Quad'
			eNodMat_ = fscanf(fid, '%d %d %d %d', [4, numEles_])'; 
		case 'Tri'
			eNodMat_ = fscanf(fid, '%d %d %d', [3, numEles_])'; 
	end

	startReadingLoads = fscanf(fid, '%s %s', 2); 
	if ~strcmp(startReadingLoads, 'NodeForces:'), error('Un-supported Data!'); end
	numLoadedNodes = fscanf(fid, '%d', 1);
	if numLoadedNodes>0, loadingCond_ = fscanf(fid, '%d %e %e', [3, numLoadedNodes])'; else, loadingCond_ = []; end
    
	startReadingFixations = fscanf(fid, '%s %s', 2);
    if ~strcmp(startReadingFixations, 'FixedNodes:'), error('Un-supported Data!'); end
	numFixedNodes = fscanf(fid, '%d', 1);
	if numFixedNodes>0, fixingCond_ = fscanf(fid, '%d', [1, numFixedNodes])'; else, fixingCond_ = []; end
    
	startReadingStress = fscanf(fid, '%s %s', 2); 
	if ~strcmp(startReadingStress, 'CartesianStress:'), error('Un-supported Data!'); end
	numStressComponents = fscanf(fid, '%d', 1);
	cartesianStressField_ = fscanf(fid, '%e %e %e', [3, numStressComponents])';		
	fclose(fid);
		
	%%Extract Boundary Element Info.
	switch meshType_
		case 'Quad'
			numNodesPerEle_ = 4;
			numEdgesPerEle_ = 4;
			eleEdges = [1 2; 2 3; 3 4; 4 1];
			edgeDirOrder = [2 3 4 1];
		case 'Tri'
			numNodesPerEle_ = 3;
			numEdgesPerEle_ = 3;
			eleEdges = [1 2; 2 3; 3 1];
			edgeDirOrder = [2 3 1];
	end
	
	boundingBox_ = [min(nodeCoords_, [], 1); max(nodeCoords_, [], 1)];	
	[boundaryEdgeNodMat_, boundaryNodes_, nodState_] = ExtractBoundaryInformation();
	%%Extract Silhouette for Vis.
	silhouetteStruct_.vertices = nodeCoords_;
	silhouetteStruct_.faces = boundaryEdgeNodMat_;

	%%element centroids
	eleNodCoordListX = nodeCoords_(:,1); eleNodCoordListX = eleNodCoordListX(eNodMat_);
	eleNodCoordListY = nodeCoords_(:,2); eleNodCoordListY = eleNodCoordListY(eNodMat_);
	eleCentroidList_ = [sum(eleNodCoordListX,2) sum(eleNodCoordListY,2)]/numNodesPerEle_;	
	
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
	
	iEleStruct = struct('edgeCentres', [], 'edgeNormals', []); %%pure-Quad
	eleStruct_ = repmat(iEleStruct, numEles_, 1);
	for ii=1:numEles_
		iNodes = eNodMat_(ii,:);
		iEleVertices = nodeCoords_(iNodes, :);
		edgeDirVec = iEleVertices(edgeDirOrder,:) - iEleVertices;
		edgeDirVec = edgeDirVec ./ vecnorm(edgeDirVec,2,2);
		aveNormal = [-edgeDirVec(:,2) edgeDirVec(:,1)];			
		tmp = iEleStruct;			
		%% tmp.edgeNormals = aveNormal;
		%% in case the node orderings on each element face are not constant
		tmp.edgeCentres = (iEleVertices + iEleVertices(edgeDirOrder,:))/2;
		iEleCt = eleCentroidList_(ii,:);
		refVecs = iEleCt - tmp.edgeCentres; refVecs = refVecs ./ vecnorm(refVecs,2,2);
		dirEval = acos(sum(refVecs .* aveNormal, 2));
		dirDes = ones(numEdgesPerEle_,1); dirDes(dirEval<pi/2) = -1;
		edgeNormals = dirDes .* aveNormal;
		tmp.edgeNormals = edgeNormals;
		eleStruct_(ii) = tmp;
	end
	
	%% Evaluate Element Sizes
	tmpSizeList = zeros(numEdgesPerEle_,numEles_);
	for ii=1:numEles_
		tmpSizeList(:,ii) = vecnorm(eleCentroidList_(ii,:)-eleStruct_(ii).edgeCentres,2,2);
	end
	eleSizeList_ = 2*min(tmpSizeList,[],1)';	
end

function [bEdgeNodMat, boundaryNodes, nodState] = ExtractBoundaryInformation()
	global meshType_;
    global numNodes_;
	global numEles_;
	global eNodMat_;
	global numEdgesPerEle_;
	switch meshType_
		case 'Quad'
			edgeIndices = eNodMat_(:, [1 2  2 3  3 4  4 1])';		
		case 'Tri'
			edgeIndices = eNodMat_(:, [1 2  2 3  3 1])';					
	end
	edgeIndices = reshape(edgeIndices(:), 2, numEdgesPerEle_*numEles_)';	
	tmp = sort(edgeIndices,2);
	[uniqueEdges, ia, ~] = unique(tmp, 'stable', 'rows');
	leftEdgeIDs = (1:numEdgesPerEle_*numEles_)'; leftEdgeIDs = setdiff(leftEdgeIDs, ia);
	leftEdges = tmp(leftEdgeIDs,:);
	[~, boundaryEdgesIDsInUniqueEdges] = setdiff(uniqueEdges, leftEdges, 'rows');
	bEdgeNodMat = edgeIndices(ia(boundaryEdgesIDsInUniqueEdges),:);
	boundaryNodes = unique(bEdgeNodMat);
	nodState = zeros(numNodes_,1);
	nodState(boundaryNodes) = 1;	
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Topology Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TensorFieldTopologyAnalysis()
	global potentialDegenerateElements_;
	global degeneratePoints_;
	global majorPSLpool_;
	global minorPSLpool_;
	
	%%0. Prepare
	PrepareTA();
	
	%%1. identify potential 'degenerate' elements
	ExtractPotentialDegenerateElements();

	%%2. identify degenerate points
	IdentifyingDegeneratePoints(potentialDegenerateElements_);	
	
	%%3. post-processing degenerate points
	degeneratePoints_ = PostProcessDegeneratePoints();
	
	%%4. get the topology skeletons
	ComputeTopologicalSkeletons();	
	
	%%5. write into PSL pools
	if length(degeneratePoints_)
		for ii=1:length(degeneratePoints_)
			if 0==degeneratePoints_(ii).eleIndex, continue; end
			for jj=1:length(degeneratePoints_(ii).majorSkeletons)
				majorPSLpool_(end+1,1) = degeneratePoints_(ii).majorSkeletons(jj);
			end
			for jj=1:length(degeneratePoints_(ii).minorSkeletons)
				minorPSLpool_(end+1,1) = degeneratePoints_(ii).minorSkeletons(jj);
			end
		end
	end	
end

function [detJ, invJ] = CalcJacobi()
	%% only for 1st-order quad element
	global meshType_;
	global eNodMat_;
	global numEles_;
	global nodeCoords_;
	global deShapeFuncs_;
	
	nEND = 2;
	switch meshType_
		case 'Quad'
			nEGIP = 4;
		case 'Tri'
			nEGIP = 3;
	end
	
	iInvJ = struct('arr', sparse(nEND*nEGIP,nEND*nEGIP));
	iDetJ = zeros(nEGIP,1);
	invJ = repmat(iInvJ, numEles_, 1);
	detJ = repmat(iDetJ,1,numEles_);
	for ii=1:numEles_
		probeEleNods = nodeCoords_(eNodMat_(ii,:)',:);
		for kk=1:nEGIP
			Jac = deShapeFuncs_(nEND*(kk-1)+1:nEND*kk,:)*probeEleNods;
			iInvJ.arr(nEND*(kk-1)+1:nEND*kk, nEND*(kk-1)+1:nEND*kk) = inv(Jac);
			iDetJ(kk) = det(Jac);	
		end
		invJ(ii) = iInvJ;
		detJ(:,ii) = iDetJ;
	end	
end

function ComputeTopologicalSkeletons()
	global nodeCoords_; 
	global eNodMat_;
    global boundaryElements_;
	global cartesianStressField_;
	global degeneratePoints_;
	global permittedMaxAdjacentTangentAngleDeviation_;
	global invJ_;
	
	%%1. get the derivatives of cartesian stresses at degenerate points with respect to the cartesian coordinates
	%%  compute rotational invariant
	permittedMaxAdjacentTangentAngleDeviation_ = 6;
	for ii=1:length(degeneratePoints_)
		s = degeneratePoints_(ii).paraCoord(1); t = degeneratePoints_(ii).paraCoord(2);
		dShape = DeShapeFunction(s, t);
		eleCoords = nodeCoords_(eNodMat_(degeneratePoints_(ii).eleIndex,:)',:);
		eleNodeCartesionStresses = cartesianStressField_(eNodMat_(degeneratePoints_(ii).eleIndex,:)',:);
		dN2dPhyC = invJ_(degeneratePoints_(ii).eleIndex).arr(1:2,1:2)*dShape;	
		degeneratePoints_(ii).stress2phy = [(dN2dPhyC(1,:)*eleNodeCartesionStresses )' (dN2dPhyC(2,:)*eleNodeCartesionStresses )' ];
		
		a = (degeneratePoints_(ii).stress2phy(1,1) - degeneratePoints_(ii).stress2phy(2,1))/2; 
		b = (degeneratePoints_(ii).stress2phy(1,2) - degeneratePoints_(ii).stress2phy(2,2))/2; 
		c = degeneratePoints_(ii).stress2phy(3,1); 
		d = degeneratePoints_(ii).stress2phy(3,2);
		degeneratePoints_(ii).abcd = [a b c d];
		degeneratePoints_(ii).delta = a*d - b*c; %% Negative -> Trisector; Positive -> Wedge
		
		rawRoots = roots([d (c+2*b) (2*a-d) -c]);
		degeneratePoints_(ii).tangentList = rawRoots(0==imag(rawRoots));		
	end
	
	%%1.1 Exclude Degenerate Points located in the Boundary Elements
	%% "An Experience-based Solution to Ease up the Numerical Instability Caused by the POTENTIAL  Jaggy Boundary of Cartesian Mesh"
	elementsWithDegeneratePoints = [degeneratePoints_.eleIndex];
	[~, boundaryElementIndicesWithDegeneratePoints] = intersect(elementsWithDegeneratePoints, boundaryElements_);
	degeneratePoints_(boundaryElementIndicesWithDegeneratePoints) = [];
	
	%%2. get the topology skeletons
	for ii=1:length(degeneratePoints_)
		degeneratePoints_(ii).majorSkeletons = repmat(degeneratePoints_(ii).majorSkeletons, length(degeneratePoints_(ii).tangentList), 1);
		degeneratePoints_(ii).minorSkeletons = repmat(degeneratePoints_(ii).minorSkeletons, length(degeneratePoints_(ii).tangentList), 1);
		for jj=1:length(degeneratePoints_(ii).tangentList)
			seed = [degeneratePoints_(ii).eleIndex degeneratePoints_(ii).phyCoord];
			iniDir = [1 degeneratePoints_(ii).tangentList(jj) ]; iniDir = iniDir/norm(iniDir);
			% [majorPSL, minorPSL] = GeneratePrincipalStressLines(seed, iniDir, []);
			majorPSL = Have1morePSL2D(seed, 'MAJOR', iniDir);
			minorPSL = Have1morePSL2D(seed, 'MINOR', iniDir);
			dis0 = degeneratePoints_(ii).phyCoord - majorPSL.phyCoordList(1,:); dis0 = norm(dis0);
			dis1 = degeneratePoints_(ii).phyCoord - majorPSL.phyCoordList(end,:); dis1 = norm(dis1);		
			if dis0 < dis1
				majorPSL.eleIndexList = majorPSL.eleIndexList(majorPSL.midPointPosition:majorPSL.length,:);
				majorPSL.phyCoordList = majorPSL.phyCoordList(majorPSL.midPointPosition:majorPSL.length,:);
				majorPSL.principalStressList = majorPSL.principalStressList(majorPSL.midPointPosition:majorPSL.length,:);
				majorPSL.midPointPosition = 1;
				majorPSL.length = length(majorPSL.eleIndexList);
			else
				majorPSL.eleIndexList = flip(majorPSL.eleIndexList(1:majorPSL.midPointPosition,:));
				majorPSL.phyCoordList = flip(majorPSL.phyCoordList(1:majorPSL.midPointPosition,:));
				majorPSL.principalStressList = flip(majorPSL.principalStressList(1:majorPSL.midPointPosition,:)); 
				majorPSL.midPointPosition = 1;
				majorPSL.length = length(majorPSL.eleIndexList);
			end
			degeneratePoints_(ii).majorSkeletons(jj) = majorPSL;
			
			dis0 = degeneratePoints_(ii).phyCoord - minorPSL.phyCoordList(1,:); dis0 = norm(dis0);
			dis1 = degeneratePoints_(ii).phyCoord - minorPSL.phyCoordList(end,:); dis1 = norm(dis1);
			if dis0 < dis1
				minorPSL.eleIndexList = minorPSL.eleIndexList(minorPSL.midPointPosition:minorPSL.length,:);
				minorPSL.phyCoordList = minorPSL.phyCoordList(minorPSL.midPointPosition:minorPSL.length,:);
				minorPSL.principalStressList = minorPSL.principalStressList(minorPSL.midPointPosition:minorPSL.length,:);
				minorPSL.midPointPosition = 1; minorPSL.length = length(minorPSL.eleIndexList);
			else
				minorPSL.eleIndexList = flip(minorPSL.eleIndexList(1:minorPSL.midPointPosition,:));
				minorPSL.phyCoordList = flip(minorPSL.phyCoordList(1:minorPSL.midPointPosition,:));
				minorPSL.principalStressList = flip(minorPSL.principalStressList(1:minorPSL.midPointPosition,:)); 
				minorPSL.midPointPosition = 1; minorPSL.length = length(minorPSL.eleIndexList);
			end		
			degeneratePoints_(ii).minorSkeletons(jj) = minorPSL;	
		end
	end
end

function IdentifyingDegeneratePoints(potentialElements)	
	global consideredDegenerateElements_; 
	global numConsideredDegenerateElements_; 
	global eNodMat_;
	global cartesianStressField_;
	global paraCoordListDegeEles_;
	
	consideredDegenerateElements_ = potentialElements(:);
	numConsideredDegenerateElements_ = length(consideredDegenerateElements_);	
	paraCoordListDegeEles_ = zeros(numConsideredDegenerateElements_,2);
	for ii=1:numConsideredDegenerateElements_
		iEleStress = cartesianStressField_(eNodMat_(potentialElements(ii),:)',:); 
		v1 = DiscriminantConstraintFuncs(iEleStress);
		[paraCoord, ~, ~, ~] = NewtonIteration(v1, zeros(1,size(v1,2)));
		paraCoordListDegeEles_(ii,:) = paraCoord;
	end
end

function extractedDegeneratePoints = PostProcessDegeneratePoints()
	global consideredDegenerateElements_;
	global paraCoordListDegeEles_;
	global numConsideredDegenerateElements_;
	global eNodMat_;
	global nodeCoords_;
	global cartesianStressField_;
	global thresholdDPE_;   
	global eleSizeList_;
	extractedDegeneratePoints = DegeneratePointStruct();
	RF = 1.0e-6; %%relaxation factor
	phyCoordList = [];
	index = 0;
	for ii=1:numConsideredDegenerateElements_	
		paraCoord = paraCoordListDegeEles_(ii,:);
		if abs(paraCoord(1))>RF+1 || abs(paraCoord(2))>RF+1, continue; end
		iDegePot = DegeneratePointStruct();
		iDegePot.eleIndex = consideredDegenerateElements_(ii);
		iDegePot.paraCoord = paraCoord;
		tarEleNodeIndices = eNodMat_(iDegePot.eleIndex,:)';
		iNodeCoord = nodeCoords_(tarEleNodeIndices,:);
		shapeFuncs = ShapeFunction(iDegePot.paraCoord(1), iDegePot.paraCoord(2));
		iDegePot.phyCoord = shapeFuncs*iNodeCoord;
		phyCoordList(end+1,:) = iDegePot.phyCoord;	
		iEleStress = cartesianStressField_(tarEleNodeIndices,:);
		iDegePot.cartesianStress = shapeFuncs*iEleStress;
		ps = ComputePrincipalStress2D(iDegePot.cartesianStress);
		iDegePot.principalStress = ps([4 1]);
		directDegenerancyMetric = abs(iDegePot.principalStress(1)-iDegePot.principalStress(2)) / ...
				abs(iDegePot.principalStress(1)+iDegePot.principalStress(2));	
		if directDegenerancyMetric>thresholdDPE_, continue; end
		iDegePot.directDegenerancyExtentMetric = directDegenerancyMetric;
		index = index + 1;	
		extractedDegeneratePoints(index,1) = iDegePot;				
	end
	if 0 == extractedDegeneratePoints(1).eleIndex, extractedDegeneratePoints(1) = []; return; end %% There is no degenerate point
	
	%%Merge the closest degenerate points
	if numel(extractedDegeneratePoints)>1
		degeneratePoints2Bmerged = [];
		degeEleList = extractedDegeneratePoints(1).eleIndex;
		degeElePosList = extractedDegeneratePoints(1).phyCoord;
		for ii=2:numel(extractedDegeneratePoints)
			iDegeEle = extractedDegeneratePoints(ii).eleIndex;
			iDegeElePos = extractedDegeneratePoints(ii).phyCoord;
			iMetric1 = isempty(setdiff(iDegeEle, degeEleList));
			[minVal, minValPos] = min(vecnorm(iDegeElePos-degeElePosList,2,2));
			iMetric2 = minVal < (eleSizeList_(degeEleList(minValPos)) + eleSizeList_(iDegeEle))/2/2;
			if iMetric1 || iMetric2
				degeneratePoints2Bmerged(end+1,1) = ii;
			else
				degeEleList(end+1,1) = iDegeEle;
				degeElePosList(end+1,:) = iDegeElePos;
			end
		end	
		extractedDegeneratePoints(degeneratePoints2Bmerged) = [];
	end
end

function [paraCoordinates, res, opt, index] = NewtonIteration(vtxVec, target)
	%% solving a nonlinear system of equasions by Newton-Rhapson's method
	%%	f1(s,t) = tar1
	%%	f2(s,t) = tar2
	opt = 0;
	normTar = norm(target);
	errThreshold = 1.0e-10; RF = 100*errThreshold;	
	s = -0.0; t = -0.0; maxIts = 150;
	index = 0;
	
	for ii=1:maxIts
		index = index+1;
		c0 = ShapeFunction(s, t)';
		dShape = DeShapeFunction(s, t);
		dns = dShape(1,:)';
		dnt = dShape(2,:)';
		d2Shape = De2ShapeFunction(s, t);
		dnss = d2Shape(1,:)';
		dntt = d2Shape(2,:)';
		dnst = d2Shape(3,:)';
		
		q = vtxVec' * c0;
		dqs = vtxVec' * dns;
		dqt = vtxVec' * dnt;

		dfdv1 = [dqs';dqt'];
		b = dfdv1*(q-target');
		if 0==normTar
			res = norm(q-target');
		else
			res = norm(b);
		end
		if res < errThreshold, break; end			
		
		dfdss = vtxVec'*dnss;
		dfdtt = vtxVec'*dntt;
		dfdst = vtxVec'*dnst;
		A11 = dfdss' * (q-target') + norm(dqs)^2;
		A22 = dfdtt' * (q-target') + norm(dqt)^2;
		A12 = dfdst' * (q-target') + dqs'*dqt;
		A21 = A12;
		A = [A11 A12; A21 A22]; x = A\(-b);		
		s = s + x(1); t = t + x(2);	
	end
	if res <= errThreshold && abs(s)<=RF+1 && abs(t)<=RF+1
		opt = 1;
	end
	paraCoordinates = [s t];
end

function val = DegeneratePointStruct()
	vec = struct('ith', 0, 'length', 0, 'vec', [], 'index', []);
	PSL = PrincipalStressLineStruct();
	val = struct(	...
		'eleIndex',							0,	...
		'paraCoord',						[],	...		
		'phyCoord',							[],	...
		'cartesianStress',					[], ...
		'principalStress',					[],	...
		'directDegenerancyExtentMetric', 	[], ...
		'tangentList',						[],	...
		'stress2phy',						[],	...
		'abcd',								[],	...
		'delta',							0,	...
		'majorSkeletons',					PSL,...
		'minorSkeletons',					PSL,...
		'separatrices',						vec	...
	);
end

function PrepareTA()
	global detJ_;
	global invJ_;	
	global deShapeFuncs_;
	
	[s, t, ~] = GaussianIntegral();
	deShapeFuncs_ = DeShapeFunction(s', t');
	[detJ_, invJ_] = CalcJacobi();	
end

function ExtractPotentialDegenerateElements()
	global eNodMat_;
	global numEles_;
	global cartesianStressField_;
	global potentialDegenerateElements_;
    global excludeDegeneratePointsOnBoundary_;
	potentialDegenerateElements_ = [];
	for ii=1:numEles_
		eleStress = cartesianStressField_(eNodMat_(ii,:)',:);
		opt = DegenrationMeasure(eleStress);
		if 1==opt, potentialDegenerateElements_(end+1,1) = ii; end			
    end

	%%Exclude the Potential Degenerate Points near the Boundary, which might be caused by the jaggy Cartesian Mesh
    if excludeDegeneratePointsOnBoundary_
        boundaryElements = ExtractBoundaryElements(excludeDegeneratePointsOnBoundary_);
        [~, compactPEWDP] = setdiff(potentialDegenerateElements_, boundaryElements);
        potentialDegenerateElements_ = potentialDegenerateElements_(compactPEWDP);        
    end    
end


function eles2Bexcluded = ExtractBoundaryElements(numLayerEles2Bexcluded)
	global nodStruct_; 
	global boundaryElements_;
	global eNodMat_;
	
	numLayerEles2Bexcluded = round(numLayerEles2Bexcluded);
	if numLayerEles2Bexcluded<1
		eles2Bexcluded = [];
	else
		eles2Bexcluded = boundaryElements_;
		idx = 2;
		while idx<=numLayerEles2Bexcluded
			iAllNodes = eNodMat_(eles2Bexcluded,:); iAllNodes = iAllNodes(:);
			iAllEles = nodStruct_(iAllNodes);
			iAllEles = [iAllEles.adjacentEles];
			iAllEles = unique(iAllEles);
			eles2Bexcluded = iAllEles(:);
		end
	end
end

function opt = DegenrationMeasure(tar)
	global meshType_;
	discriminants = DiscriminantConstraintFuncs(tar);
	v1 = discriminants(:,1);
	v2 = discriminants(:,2);
	
	switch meshType_
		case 'Quad'
			bool1_1 = v1(1)>0 && v1(2)>0 && v1(3)>0 && v1(4)>0;
			bool1_2 = v1(1)<0 && v1(2)<0 && v1(3)<0 && v1(4)<0;	

			bool2_1 = v2(1)>0 && v2(2)>0 && v2(3)>0 && v2(4)>0;
			bool2_2 = v2(1)<0 && v2(2)<0 && v2(3)<0 && v2(4)<0;			
		case 'Tri'
			bool1_1 = v1(1)>0 && v1(2)>0 && v1(3)>0;
			bool1_2 = v1(1)<0 && v1(2)<0 && v1(3)<0;	

			bool2_1 = v2(1)>0 && v2(2)>0 && v2(3)>0;
			bool2_2 = v2(1)<0 && v2(2)<0 && v2(3)<0;			
	end

	bool1 = bool1_1 || bool1_2;
	bool2 = bool2_1 || bool2_2;
	
	if bool1 || bool2, opt = 0; else, opt = 1; end	
end

function discriminants = DiscriminantConstraintFuncs(eleStress)
	discriminants = [eleStress(:,1)-eleStress(:,2), eleStress(:,3)];
end

function N = ShapeFunction(s, t)
	%% 1st-order Quadrilateral plane element (4 nodes and 4 Gaussian integral points)
	%				   	   __s (parametric coordinate system)
	%				  	 /-t
	%				*4			*3
	%			*1			*2
	%
	%				nodes

	%% 1st-order Triangular plane element (3 nodes and 3 Gaussian integral points)
	%			3						*3
	%			|  \			 	
	%	y   	|    \						t
	%	|__x	|	 	\ 					|__s
	%			|	       \	 	
	%			1-----------2			*1			*2
	%			Node Ordering		Gauss IPs Ordering
	global meshType_;
	s = s(:);
	t = t(:);
	switch meshType_
		case 'Quad'
			N = zeros(size(s,1), 4);
			N(:,1) = 0.25*(1-s).*(1-t);
			N(:,2) = 0.25*(1+s).*(1-t);
			N(:,3) = 0.25*(1+s).*(1+t);
			N(:,4) = 0.25*(1-s).*(1+t);			
		case 'Tri'
			N(:,1) = 1-s-t;
			N(:,2) = s;
			N(:,3) = t;		
	end

end

function dN = DeShapeFunction(s, t)
	global meshType_;
	s = s(:);
	t = t(:);
	switch meshType_
		case 'Quad'
			dN1ds = -(1-t); dN2ds = 1-t; 	dN3ds = 1+t; dN4ds = -(1+t);
			dN1dt = -(1-s); dN2dt = -(1+s); dN3dt = 1+s; dN4dt = 1-s;
			
			dN = zeros(2*length(s), 4);
			dN(1:2:end,:) = 0.25*[dN1ds dN2ds dN3ds dN4ds];
			dN(2:2:end,:) = 0.25*[dN1dt dN2dt dN3dt dN4dt];			
		case 'Tri'
			tmp = ones(numel(s),1);
			dN1ds = -1*tmp; dN2ds = 1*tmp; dN3ds = 0*tmp; 
			dN1dt = -1*tmp; dN2dt = 0*tmp; dN3dt = 1*tmp;
			dN = zeros(2*length(s), 3);
			dN(1:2:end,:) = [dN1ds dN2ds dN3ds];
			dN(2:2:end,:) = [dN1dt dN2dt dN3dt];				
	end
end

function d2Shape = De2ShapeFunction(s, t)	
	global meshType_;
	s = s(:);
	t = t(:);
	numCoord = length(s);
	switch meshType_
		case 'Quad'
			dN1dss = 0; dN2dss = 0; dN3dss = 0; dN4dss = 0;
			dN1dtt = 0; dN2dtt = 0; dN3dtt = 0; dN4dtt = 0;
			dN1dst = 0.25; dN2dst = -0.25; dN3dst = 0.25; dN4dst = -0.25;	
			
			d2Shape = repmat(zeros(size(s)),3,4);
			d2Shape(1:3:end,:) = repmat([dN1dss	dN2dss	dN3dss	dN4dss], numCoord, 1);
			d2Shape(2:3:end,:) = repmat([dN1dtt	dN2dtt	dN3dtt	dN4dtt], numCoord, 1);
			d2Shape(3:3:end,:) = repmat([dN1dst	dN2dst	dN3dst	dN4dst], numCoord, 1);		
		case 'Tri'
			dN1dss = 0; dN2dss = 0; dN3dss = 0;
			dN1dtt = 0; dN2dtt = 0; dN3dtt = 0;
			dN1dst = 0; dN2dst = 0; dN3dst = 0;
			d2Shape = repmat(s,3,3);
			d2Shape(1:3:end,:) = repmat([dN1dss	dN2dss	dN3dss], numCoord, 1);
			d2Shape(2:3:end,:) = repmat([dN1dtt	dN2dtt	dN3dtt], numCoord, 1);
			d2Shape(3:3:end,:) = repmat([dN1dst	dN2dst	dN3dst], numCoord, 1);					
	end

end

function [s, t, w] = GaussianIntegral()
	%		*4			*3
	%
	%			t
	%			|__s
	%		
	%		*1			*2
		
	%		*3
	% 	
	%			t
	%			|__s
	% 	
	%		*1			*2
	%		Gaussian point
	global meshType_;
	switch meshType_
		case 'Quad'
			s = [-0.577350269189626 0.577350269189626 0.577350269189626 -0.577350269189626]';
			t = [-0.577350269189626 -0.577350269189626 0.577350269189626 0.577350269189626]';
			w = [1.0, 1.0, 1.0, 1.0]';		
		case 'Tri'
			s = [1/6 2/3 1/6]';
			t = [1/6 1/6 2/3]';
			w = [1/6 1/6 1/6]';			
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowPSLs()
	global silhouetteStruct_;
	global majorPSLpool_;
	global minorPSLpool_;
	global degeneratePoints_;
	
	miniLength2Bshown = 5;
	for jj=1:length(majorPSLpool_)
		if majorPSLpool_(jj).length>miniLength2Bshown
			plot(majorPSLpool_(jj).phyCoordList(:,1), majorPSLpool_(jj).phyCoordList(:,2), ...
				'-', 'color', [252 141 98]/255, 'LineWidth', 3); hold('on');		
		end
	end	
	for jj=1:length(minorPSLpool_)
		if minorPSLpool_(jj).length>miniLength2Bshown
			plot(minorPSLpool_(jj).phyCoordList(:,1), minorPSLpool_(jj).phyCoordList(:,2), ...
				'-', 'color', [102 194 165]/255, 'LineWidth', 3); hold('on');
		end
	end	

	if length(degeneratePoints_) > 0
		for ii=1:length(degeneratePoints_)		
			hold('on');
            if degeneratePoints_(ii).delta > 0
				plot(degeneratePoints_(ii).phyCoord(1), degeneratePoints_(ii).phyCoord(2), 'sm', 'LineWidth', 3, 'MarkerSize', 15);
			else
				plot(degeneratePoints_(ii).phyCoord(1), degeneratePoints_(ii).phyCoord(2), 'om', 'LineWidth', 3, 'MarkerSize', 15);
			end		
		end		
	end

	%%Show silhouette
	hSilo = patch(silhouetteStruct_); hold('on');
	set(hSilo, 'FaceColor', 'None', 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 2);
	axis('equal'); axis('tight'); axis('off');	
end

function ShowProblemDescription()
	global nodeCoords_;
	global eNodMat_;
	global loadingCond_;
	global fixingCond_;
	global boundingBox_;
	
	tmp.vertices = nodeCoords_;
	tmp.faces = eNodMat_;
	
	hd = patch(tmp); hold('on');
	set(hd, 'FaceColor', [65 174 118]/255, 'EdgeColor', 'k');
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
		hold('on'); quiver(coordLoadedNodes(:,1), coordLoadedNodes(:,2), amplitudesF(:,1), ...
			amplitudesF(:,2), 0, 'Color', [255 127 0.0]/255, 'LineWidth', 2, 'MaxHeadSize', 1, 'MaxHeadSize', 1);
	end
    if ~isempty(fixingCond_)
		tarNodeCoord = nodeCoords_(fixingCond_(:,1),:);
		hold('on'); hd1 = plot(tarNodeCoord(:,1), tarNodeCoord(:,2), 'x', ...
			'color', [153 153 153]/255, 'LineWidth', 3, 'MarkerSize', 15);		
	end
	axis('equal'); axis('tight'); axis('off');	
end