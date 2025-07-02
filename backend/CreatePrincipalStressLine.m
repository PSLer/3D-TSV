function iPSL = CreatePrincipalStressLine(startPoint, tracingType)
	global tracingFuncHandle_;
	global sPoint_;
	global integrationStepLimit_;
	sPoint_ = startPoint(end-2:end);
	iPSL = PrincipalStressLineStruct();
	switch tracingType
		case 'MAJOR', psDir = [10 11 12];
		case 'MEDIUM', psDir = [6 7 8];		
		case 'MINOR', psDir = [2 3 4];
	end
	%%1. prepare for tracing			
	[eleIndex, cartesianStress, vonMisesStress, principalStress, opt] = PreparingForTracing(startPoint);
	if 0==opt, return; end
	
	%%2. tracing PSL
	startPoint = startPoint(end-2:end);
	PSLphyCoordList = startPoint;
	PSLcartesianStressList = cartesianStress;
	PSLeleIndexList = eleIndex;
	PSLvonMisesStressList = vonMisesStress;
	PSLprincipalStressList = principalStress;			
	%%2.1 along first direction (v1)		
	[phyCoordList, cartesianStressList, eleIndexList, ~, vonMisesStressList, principalStressList] = ...
		tracingFuncHandle_(startPoint, principalStress(1,psDir), eleIndex, psDir, integrationStepLimit_);		
	PSLphyCoordList = [PSLphyCoordList; phyCoordList];
	PSLcartesianStressList = [PSLcartesianStressList; cartesianStressList];
	PSLeleIndexList = [PSLeleIndexList; eleIndexList];
	PSLvonMisesStressList = [PSLvonMisesStressList; vonMisesStressList];
	PSLprincipalStressList = [PSLprincipalStressList; principalStressList];
	if size(phyCoordList,1)>0, sPoint_ = phyCoordList(end,:); end
	%%2.2 along second direction (-v1)	
	[phyCoordList, cartesianStressList, eleIndexList, ~, vonMisesStressList, principalStressList] = ...
		tracingFuncHandle_(startPoint, -principalStress(1,psDir), eleIndex, psDir, integrationStepLimit_);		
	if size(phyCoordList,1) > 1
		phyCoordList = flip(phyCoordList);
		cartesianStressList = flip(cartesianStressList);
		eleIndexList = flip(eleIndexList);
		vonMisesStressList = flip(vonMisesStressList);
		principalStressList = flip(principalStressList);
	end						
	PSLphyCoordList = [phyCoordList; PSLphyCoordList];
	PSLcartesianStressList = [cartesianStressList; PSLcartesianStressList];
	PSLeleIndexList = [eleIndexList; PSLeleIndexList];
	PSLvonMisesStressList = [vonMisesStressList; PSLvonMisesStressList];
	PSLprincipalStressList = [principalStressList; PSLprincipalStressList];
	%%2.3 finish Tracing the current major PSL	
	iPSL.midPointPosition = size(phyCoordList,1)+1;
	iPSL.length = size(PSLphyCoordList,1);
	iPSL.eleIndexList = PSLeleIndexList;
	iPSL.phyCoordList = PSLphyCoordList;
	iPSL.cartesianStressList = PSLcartesianStressList;	
	iPSL.vonMisesStressList = PSLvonMisesStressList;
	iPSL.principalStressList = PSLprincipalStressList;	
end

function [eleIndex, cartesianStress, vonMisesStress, principalStress, opt] = PreparingForTracing(startPoint)
	global nodeCoords_; 
	global eNodMat_; 
	global cartesianStressField_;
	global eleCentroidList_;
	global meshType_;

	eleIndex = 0;
	cartesianStress = 0;
	vonMisesStress = 0; 
	principalStress = 0;	
	if strcmp(meshType_, 'CARTESIAN_GRID')	
		[targetEleIndex, paraCoordinates, opt] = SearchNextIntegratingPointOnCartesianMesh(startPoint(end-2:end));	
		if 0==opt, return; end
		eleIndex = double(targetEleIndex);
		NIdx = eNodMat_(eleIndex,:)';
		eleCartesianStress = cartesianStressField_(NIdx,:);				
		cartesianStress = ElementInterpolationTrilinear(eleCartesianStress, paraCoordinates);
	else
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
	end
	vonMisesStress = ComputeVonMisesStress(cartesianStress);
	principalStress = ComputePrincipalStress(cartesianStress);		
end