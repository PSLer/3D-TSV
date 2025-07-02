function [ribbonVertices, facePatches, outlinePatches, faceColors] = ExpandPSLs2RibbonsSim(PSLs, lw, psDir, colorSrc, smoothingOpt)
	%%			RIBBON
	%%	===========================
	%%		   dir2 |	
	%%				 ---dir1
	%%			   / dir3
	%%	===========================
	%%
	ribbonVertices = [];
	facePatches = [];
	outlinePatches = [];
	faceColors = [];
	
	if isempty(PSLs), return; end
	twistThreshold = 3.5/180*pi;
	numPSLs = length(PSLs);
	
	for ii=1:numPSLs
		%%1. ribbon boundary nodes
		iPSLength = PSLs(ii).length;
		iCoordList = zeros(2*iPSLength,3);
		midPots = PSLs(ii).phyCoordList;
		
		dirVecs = PSLs(ii).principalStressList(:,psDir);
		angList = zeros(iPSLength-1,1);
		for jj=2:iPSLength
			vec0 = dirVecs(jj-1,:);
			vec1 = dirVecs(jj,:); vec2 = -vec1;
			ang1 = acos(vec0 * vec1'); ang2 = acos(vec0 * vec2');
			angList(jj-1) = ang1;
			if ang2<ang1
				angList(jj-1) = ang2;
				dirVecs(jj,:) = vec2;
			end
		end
		if smoothingOpt
			angDeviationMetric = sum(angList)/(iPSLength-1);	
			if angDeviationMetric>twistThreshold
				for jj=2:iPSLength
					vec0 = dirVecs(jj-1,:);
					vec1 = dirVecs(jj,:);
					ang1 = acos(vec0 * vec1');
					if ang1>angDeviationMetric
						dirVecs(jj,:) = vec0;
					end
				end
			end
		end				
		dirVecs = dirVecs * lw;
		
		coords1 = midPots + dirVecs;
		coords2 = midPots - dirVecs;
		iCoordList(1:2:end,:) = coords1;
		iCoordList(2:2:end,:) = coords2;
			
		%%2. create quad patches
		numExistingNodes = size(ribbonVertices,1);
		numNewlyGeneratedNodes = 2*iPSLength;
		newGeneratedNodes = numExistingNodes + (1:numNewlyGeneratedNodes);
		newGeneratedNodes = reshape(newGeneratedNodes, 2, iPSLength);
		iQuadMapFace = [newGeneratedNodes(1,1:end-1); newGeneratedNodes(2,1:end-1); ...
			newGeneratedNodes(2,2:end); newGeneratedNodes(1,2:end)];
			
		%%3. write into global ribbon info
		iFaceColorList = colorSrc(ii).arr;
		iFaceColorList = reshape(repmat(iFaceColorList, 2, 1), 2*iPSLength, 1);	
		iQuadMapOutline = [
			newGeneratedNodes(1,1:end-1) newGeneratedNodes(2,1:end-1) newGeneratedNodes(1,1) newGeneratedNodes(1,end)
			newGeneratedNodes(1,2:end)	 newGeneratedNodes(2,2:end)	  newGeneratedNodes(2,1) newGeneratedNodes(2,end)
			newGeneratedNodes(1,2:end)	 newGeneratedNodes(2,2:end)   newGeneratedNodes(2,1) newGeneratedNodes(2,end)
			newGeneratedNodes(1,1:end-1) newGeneratedNodes(2,1:end-1) newGeneratedNodes(1,1) newGeneratedNodes(1,end)
		];		
		faceColors(end+1:end+2*iPSLength,:) = iFaceColorList;
		outlinePatches(end+1:end+2*iPSLength,:) = iQuadMapOutline';
		facePatches(end+1:end+iPSLength-1,:) = iQuadMapFace';
		
		ribbonVertices(end+1:end+2*iPSLength,:) = iCoordList;
	end		
end
