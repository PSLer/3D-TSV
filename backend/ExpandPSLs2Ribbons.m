function [ribbonVertices, facePatches, faceColors, faceIndices] = ExpandPSLs2Ribbons(PSLs, lw, psDir, colorSrc)
	%%			RIBBON
	%%	===========================
	%%		   dir2 |	
	%%				 ---dir1
	%%			   / dir3
	%%	===========================
	%%
	ribbonVertices = [];
	facePatches = [];
	faceColors = [];
	faceIndices = [];
	
	if isempty(PSLs), return; end
	numPSLs = length(PSLs);
	faceIndices = struct('arr', []);
	faceIndices = repmat(faceIndices, numPSLs, 1);		
	for ii=1:numPSLs
		%%1. ribbon boundary nodes
		iPSLength = PSLs(ii).length;
		iCoordList = zeros(4*iPSLength,3);
		midPots = PSLs(ii).phyCoordList;
		
		dirVecs = PSLs(ii).principalStressList(:,psDir(1,:));	
		angList = zeros(iPSLength-1,1);
		
		dirVecsDeputy = PSLs(ii).principalStressList(:,psDir(2,:));
		angListDeputy = angList;
		for jj=2:iPSLength
			vec0 = dirVecs(jj-1,:);
			vec1 = dirVecs(jj,:); vec2 = -vec1;
			ang1 = acos(vec0 * vec1'); ang2 = acos(vec0 * vec2');
			angList(jj-1) = ang1;
			if ang2<ang1
				angList(jj-1) = ang2;
				dirVecs(jj,:) = vec2;
			end
			
			vec0Deputy = dirVecsDeputy(jj-1,:);
			vec1Deputy = dirVecsDeputy(jj,:); vec2Deputy = -vec1Deputy;
			ang1Deputy = acos(vec0Deputy * vec1Deputy'); ang2Deputy = acos(vec0Deputy * vec2Deputy');
			angListDeputy(jj-1) = ang1Deputy;
			if ang2Deputy<ang1Deputy
				angListDeputy(jj-1) = ang2Deputy;
				dirVecsDeputy(jj,:) = vec2Deputy;
			end			
		end
			
		dirVecs = dirVecs * lw;
		dirVecsDeputy = dirVecsDeputy * lw * 0.15;
		
		coordsRef1 = midPots + dirVecs;
		coordsRef2 = midPots - dirVecs;
		
		coords1 = coordsRef1 + dirVecsDeputy;
		coords2 = coordsRef2 + dirVecsDeputy;
		coords3 = coordsRef2 - dirVecsDeputy;
		coords4 = coordsRef1 - dirVecsDeputy;
		
		iCoordList(1:4:end,:) = coords1;
		iCoordList(2:4:end,:) = coords2;
		iCoordList(3:4:end,:) = coords3;
		iCoordList(4:4:end,:) = coords4;
			
		%%2. create quad patches
		numExistingNodes = size(ribbonVertices,1);
		numNewlyGeneratedNodes = 4*iPSLength;
		
		iQuadMapFace = zeros(4*(iPSLength-1)+2,4);
		tmp = 4*ones(iPSLength-1,4) .* (0:iPSLength-2)'; 
		
		%%f1
		islider = 0*(iPSLength-1);
		iQuadMapFace(1+islider:iPSLength-1+islider,:) = repmat([1 2 6 5], iPSLength-1, 1) + tmp;
		%%f2
		islider = 1*(iPSLength-1);
		iQuadMapFace(1+islider:iPSLength-1+islider,:) = repmat([4 8 7 3], iPSLength-1, 1) + tmp;
		%%f3
		islider = 2*(iPSLength-1);
		iQuadMapFace(1+islider:iPSLength-1+islider,:) = repmat([1 5 8 4], iPSLength-1, 1) + tmp;		
		%%f4
		islider = 3*(iPSLength-1);
		iQuadMapFace(1+islider:iPSLength-1+islider,:) = repmat([2 3 7 6], iPSLength-1, 1) + tmp;
		
		%%cap
		iQuadMapFace(end-1,:) = [1 4 3 2];
		iQuadMapFace(end,:) = [1 2 3 4] + 4*(iPSLength-1);
	
		newGeneratedNodes = numExistingNodes + (1:numNewlyGeneratedNodes)';
		
		
		%%3. write into global ribbon info
		iFaceColorList = colorSrc(ii).arr;
		iFaceColorList = reshape(repmat(iFaceColorList, 4, 1), 4*iPSLength, 1);	
		
		faceColors(end+1:end+4*iPSLength,:) = iFaceColorList;
		
		faceIndices(ii).arr = size(facePatches,1)+1 : size(facePatches,1)+4*(iPSLength-1)+2;
		
		facePatches(end+1:end+4*(iPSLength-1)+2,:) = newGeneratedNodes(iQuadMapFace);
		
		ribbonVertices(end+1:end+4*iPSLength,:) = iCoordList;
	end
end
