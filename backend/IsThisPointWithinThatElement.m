function opt = IsThisPointWithinThatElement(tarEleIndex, iCoord, debug)
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

	%%===============================for debug
    if debug
		global eleFaceIndices_;
		global eleSizeList_;
		if ~isempty(find(0==[refVecNorm; refVec2VerticesNorm]))
			refVec = refVec ./ refVecNorm;
			normVecs = eleStruct_(tarEleIndex).faceNormals;
			angleDevs = zeros(numFacesPerEle_,1);
			for ii=1:numFacesPerEle_
				angleDevs(ii) = acos(refVec(ii,:)*normVecs(ii,:)')/pi*180;
			end			
		end
		angleDevs
		% eleFaces = [4 3 2 1; 5 6 7 8; 1 2 6 5; 8 7 3 4; 5 8 4 1; 2 3 7 6];
		cellPatchs.faces = eleFaceIndices_;
		cellPatchs.vertices = nodeCoords_(eNodMat_(tarEleIndex,:),:);	
		% patchX = cellPatchs.vertices(:,1); patchX = patchX(cellPatchs.faces);
		% patchY = cellPatchs.vertices(:,2); patchY = patchY(cellPatchs.faces);
		% patchZ = cellPatchs.vertices(:,3); patchZ = patchZ(cellPatchs.faces);
		% patchC = zeros(size(patchZ));
		
		% featureSize = sum(vecnorm(sum(cellPatchs.vertices,1)/8 - cellPatchs.vertices,2,2))/8*2;
		featureSize = eleSizeList_(tarEleIndex);
		hdRefPot = plot3(iCoord(1), iCoord(2), iCoord(3), '*r', 'LineWidth', 3, 'MarkerSize', 15); hold on;
		% hdEle = patch(patchX',patchY',patchZ',patchC'); hold on;
		hdEle = patch(cellPatchs); hold on;
		hdCentres = plot3(iElefaceCentres(:,1), iElefaceCentres(:,2), iElefaceCentres(:,3), '+k', 'LineWidth', 3, 'MarkerSize', 15); hold on;
		normVecsToDraw = normVecs ./ vecnorm(normVecs,2,2) * featureSize;	
		hdNormals = quiver3(iElefaceCentres(:,1),iElefaceCentres(:,2),iElefaceCentres(:,3),normVecsToDraw(:,1),normVecsToDraw(:,2),normVecsToDraw(:,3)); hold on;
		refPos = repmat(iCoord,6,1);
		refVecToDraw = refVec ./ vecnorm(refVec,2,2) * featureSize/1.5;	
		hdNormals2 = quiver3(refPos(:,1),refPos(:,2),refPos(:,3),refVecToDraw(:,1),refVecToDraw(:,2),refVecToDraw(:,3)); hold on;
		set(hdEle, 'faceColor', [0.5 0.5 0.5], 'faceAlpha', 0.3);
		set(hdNormals, 'LineWidth', 3, 'Color', [255 127 0]/255);
		set(hdNormals2, 'LineWidth', 3, 'Color', [51 160 44]/255);
		xlabel('X'); ylabel('Y'); zlabel('Z');
		axis equal; axis tight; axis on;        
    end
end