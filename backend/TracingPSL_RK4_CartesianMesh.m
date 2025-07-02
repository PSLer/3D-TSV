function [phyCoordList, cartesianStressList, eleIndexList, paraCoordList, vonMisesStressList, principalStressList] = ...
			TracingPSL_RK4_CartesianMesh(startPoint, iniDir, elementIndex, typePSL, limiSteps)
	global eNodMat_;
	global cartesianStressField_;
	global tracingStepWidth_;
	global sPoint_;
	siE = 1.0e-06;
	phyCoordList = zeros(limiSteps,3);
	cartesianStressList = zeros(limiSteps,6);
	eleIndexList = zeros(limiSteps,1);
	paraCoordList = [];
	vonMisesStressList = zeros(limiSteps,1);
	principalStressList = zeros(limiSteps,12);

	%%initialize initial k1, k2, k3 and k4
	k1 = iniDir;
	midPot1 = startPoint + k1*tracingStepWidth_/2;
	index = 0;
	[elementIndex2, paraCoordinates2, bool1] = SearchNextIntegratingPointOnCartesianMesh(midPot1);
	if bool1
		%%k2
		cartesianStress2 = cartesianStressField_(eNodMat_(elementIndex2,:)', :);
		cartesianStressOnGivenPoint2 = ElementInterpolationTrilinear(cartesianStress2, paraCoordinates2);
		principalStress2 = ComputePrincipalStress(cartesianStressOnGivenPoint2);
		[k2, terminationCond] = BidirectionalFeatureProcessing(k1, principalStress2(typePSL));
		midPot2 = startPoint + tracingStepWidth_*k2/2;
		%%k3
		[elementIndex3, paraCoordinates3, bool1] = SearchNextIntegratingPointOnCartesianMesh(midPot2);
		if bool1
			cartesianStress3 = cartesianStressField_(eNodMat_(elementIndex3,:)', :);
			cartesianStressOnGivenPoint3 = ElementInterpolationTrilinear(cartesianStress3, paraCoordinates3);
			principalStress3 = ComputePrincipalStress(cartesianStressOnGivenPoint3);
			[k3, terminationCond] = BidirectionalFeatureProcessing(k1, principalStress3(typePSL));
			midPot3 = startPoint + tracingStepWidth_*k3;
			%%k4
			[elementIndex4, paraCoordinates4, bool1] = SearchNextIntegratingPointOnCartesianMesh(midPot3);
			if bool1
				cartesianStress4 = cartesianStressField_(eNodMat_(elementIndex4,:)', :);
				cartesianStressOnGivenPoint4 = ElementInterpolationTrilinear(cartesianStress4, paraCoordinates4);
				principalStress4 = ComputePrincipalStress(cartesianStressOnGivenPoint4);
				[k4, terminationCond] = BidirectionalFeatureProcessing(k1, principalStress4(typePSL));
				nextPoint = startPoint + tracingStepWidth_ * (k1 + 2*k2 + 2*k3 + k4)/6;
				[elementIndex, paraCoordinates, bool1] = SearchNextIntegratingPointOnCartesianMesh(nextPoint);
				while bool1
					index = index + 1; if index > limiSteps, index = index-1; break; end						
					cartesianStress = cartesianStressField_(eNodMat_(elementIndex,:)', :);
					cartesianStressOnGivenPoint = ElementInterpolationTrilinear(cartesianStress, paraCoordinates);
					vonMisesStress = ComputeVonMisesStress(cartesianStressOnGivenPoint);
					principalStress = ComputePrincipalStress(cartesianStressOnGivenPoint);
					evs = principalStress([1 5 9]);
					if min([abs((evs(1)-evs(2))/2/(evs(1)+evs(2))) abs((evs(3)-evs(2))/2/(evs(3)+evs(2)))])<siE, index = index-1; break; end %%degenerate point							
					%%k1
					[k1, terminationCond] = BidirectionalFeatureProcessing(iniDir, principalStress(typePSL));
					if ~terminationCond, index = index-1; break; end
					midPot1 = nextPoint + k1*tracingStepWidth_/2;
					%%k2
					[elementIndex2, paraCoordinates2, bool1] = SearchNextIntegratingPointOnCartesianMesh(midPot1);
					if ~bool1, index = index-1; break; end
					cartesianStress2 = cartesianStressField_(eNodMat_(elementIndex2,:)', :);
					cartesianStressOnGivenPoint2 = ElementInterpolationTrilinear(cartesianStress2, paraCoordinates2);
					principalStress2 = ComputePrincipalStress(cartesianStressOnGivenPoint2);
					[k2, ~] = BidirectionalFeatureProcessing(iniDir, principalStress2(typePSL));						
					midPot2 = nextPoint + k2*tracingStepWidth_/2;
					%%k3
					[elementIndex3, paraCoordinates3, bool1] = SearchNextIntegratingPointOnCartesianMesh(midPot2);
					if ~bool1, index = index-1; break; end
					cartesianStress3 = cartesianStressField_(eNodMat_(elementIndex3,:)', :);
					cartesianStressOnGivenPoint3 = ElementInterpolationTrilinear(cartesianStress3, paraCoordinates3);
					principalStress3 = ComputePrincipalStress(cartesianStressOnGivenPoint3);
					[k3, ~] = BidirectionalFeatureProcessing(iniDir, principalStress3(typePSL));						
					midPot3 = nextPoint + k3*tracingStepWidth_;
					%%k4
					[elementIndex4, paraCoordinates4, bool1] = SearchNextIntegratingPointOnCartesianMesh(midPot3);
					if ~bool1, index = index-1; break; end
					cartesianStress4 = cartesianStressField_(eNodMat_(elementIndex4,:)', :);
					cartesianStressOnGivenPoint4 = ElementInterpolationTrilinear(cartesianStress4, paraCoordinates4);
					principalStress4 = ComputePrincipalStress(cartesianStressOnGivenPoint4);
					[k4, ~] = BidirectionalFeatureProcessing(iniDir, principalStress4(typePSL));
					%%store	
					iniDir = k1;
					phyCoordList(index,:) = nextPoint;
					cartesianStressList(index,:) = cartesianStressOnGivenPoint;
					eleIndexList(index,:) = elementIndex;			
					vonMisesStressList(index,:) = vonMisesStress;
					principalStressList(index,:) = principalStress;
					%%next point
					nextPoint = nextPoint + tracingStepWidth_ * (k1 + 2*k2 + 2*k3 + k4)/6;
					if norm(sPoint_-nextPoint)<tracingStepWidth_, break; end
					[elementIndex, paraCoordinates, bool1] = SearchNextIntegratingPointOnCartesianMesh(nextPoint);						
				end
			end
		end
	end
	phyCoordList = phyCoordList(1:index,:);
	cartesianStressList = cartesianStressList(1:index,:);
	eleIndexList = eleIndexList(1:index,:);
	vonMisesStressList = vonMisesStressList(1:index,:);
	principalStressList = principalStressList(1:index,:);	
end
