function ExportResult(fileName)	
	global boundingBox_;
	global meshType_;
	global minimumEpsilon_;
	global majorPSLpool_;
	global mediumPSLpool_;
	global minorPSLpool_;
	global majorHierarchy_;
	global mediumHierarchy_;	
	global minorHierarchy_; 
	global PSLsAppearanceOrder_;
	global surfaceMeshNodeCoords_;
	global surfaceMeshElements_;

	lineWidthTube = min([1.5*minimumEpsilon_/5, min(boundingBox_(2,:)-boundingBox_(1,:))/10]);
	
	fid = fopen(fileName, 'w');
	outPutFormat1 = '%16.6e';
	outPutFormat2 = ' %.6e';
	outPutFormat3 = ' %.6f';
	outPutFormat = outPutFormat2;
	%%1. write PSLs
	iColor = struct('arr', []);
	%%1.1 major
	numMajorPSLs = length(majorPSLpool_);
	majorPSLsGlobalAppearanceOrder = find(1==PSLsAppearanceOrder_(:,1));
	fprintf(fid, '%s', '#Major'); fprintf(fid, ' %d\n', numMajorPSLs);
	for ii=1:numMajorPSLs
		iPSL = majorPSLpool_(ii);
		iSeed = iPSL.phyCoordList(iPSL.midPointPosition,:);
		pslCoords = iPSL.phyCoordList; 
		pslCoords = reshape(pslCoords', numel(pslCoords), 1)';
		iColor.arr = ones(1,iPSL.length);
		[ribbonCoordsUnsmoothed, ~, ~, ~] = ExpandPSLs2RibbonsSim(iPSL, lineWidthTube, [6 7 8], iColor, 0);
		ribbonCoordsUnsmoothed = reshape(ribbonCoordsUnsmoothed', numel(ribbonCoordsUnsmoothed), 1)';
		[ribbonCoordsSmoothed, ~, ~, ~] = ExpandPSLs2RibbonsSim(iPSL, lineWidthTube, [6 7 8], iColor, 1);
		ribbonCoordsSmoothed = reshape(ribbonCoordsSmoothed', numel(ribbonCoordsSmoothed), 1)';
		stressScalarFields4ColorCoding = [iPSL.principalStressList(:,9) iPSL.vonMisesStressList iPSL.cartesianStressList];
		
		fprintf(fid, '%d %.6f %.6f %.6f %.6f', [iPSL.length majorHierarchy_(ii,:)]);	
		fprintf(fid, ' %d %.6f %.6f %.6f\n', [majorPSLsGlobalAppearanceOrder(ii) iSeed]);	
		fprintf(fid, outPutFormat, pslCoords); fprintf(fid, '\n');
		fprintf(fid, outPutFormat, ribbonCoordsUnsmoothed); fprintf(fid, '\n');
		fprintf(fid, outPutFormat, ribbonCoordsSmoothed); fprintf(fid, '\n');
		for jj=1:8
			fprintf(fid, outPutFormat, stressScalarFields4ColorCoding(:,jj)'); fprintf(fid, '\n');
		end
	end
	%%1.2 medium
	numMediumPSLs = length(mediumPSLpool_);
	mediumPSLsGlobalAppearanceOrder = find(2==PSLsAppearanceOrder_(:,1));
	fprintf(fid, '%s', '#Medium'); fprintf(fid, ' %d\n', numMediumPSLs);
	for ii=1:numMediumPSLs
		iPSL = mediumPSLpool_(ii);
		iSeed = iPSL.phyCoordList(iPSL.midPointPosition,:);
		pslCoords = iPSL.phyCoordList; 
		pslCoords = reshape(pslCoords', numel(pslCoords), 1)';
		iColor.arr = ones(1,iPSL.length);
		[ribbonCoordsUnsmoothed, ~, ~, ~] = ExpandPSLs2RibbonsSim(iPSL, lineWidthTube, [2 3 4], iColor, 0);
		ribbonCoordsUnsmoothed = reshape(ribbonCoordsUnsmoothed', numel(ribbonCoordsUnsmoothed), 1)';
		[ribbonCoordsSmoothed, ~, ~, ~] = ExpandPSLs2RibbonsSim(iPSL, lineWidthTube, [2 3 4], iColor, 1);
		ribbonCoordsSmoothed = reshape(ribbonCoordsSmoothed', numel(ribbonCoordsSmoothed), 1)';
		stressScalarFields4ColorCoding = [iPSL.principalStressList(:,5) iPSL.vonMisesStressList iPSL.cartesianStressList];

		fprintf(fid, '%d %.6f %.6f %.6f %.6f', [iPSL.length mediumHierarchy_(ii,:)]);		
		fprintf(fid, ' %d %.6f %.6f %.6f\n', [mediumPSLsGlobalAppearanceOrder(ii) iSeed]);		
		fprintf(fid, outPutFormat, pslCoords); fprintf(fid, '\n');
		fprintf(fid, outPutFormat, ribbonCoordsUnsmoothed); fprintf(fid, '\n');
		fprintf(fid, outPutFormat, ribbonCoordsSmoothed); fprintf(fid, '\n');
		for jj=1:8
			fprintf(fid, outPutFormat, stressScalarFields4ColorCoding(:,jj)'); fprintf(fid, '\n');
		end
	end		
	%%1.3 minor
	numMinorPSLs = length(minorPSLpool_);
	minorPSLsGlobalAppearanceOrder = find(3==PSLsAppearanceOrder_(:,1));
	fprintf(fid, '%s', '#Minor'); fprintf(fid, ' %d\n', numMinorPSLs);
	for ii=1:numMinorPSLs
		iPSL = minorPSLpool_(ii);
		iSeed = iPSL.phyCoordList(iPSL.midPointPosition,:);
		pslCoords = iPSL.phyCoordList; 
		pslCoords = reshape(pslCoords', numel(pslCoords), 1)';
		iColor.arr = ones(1,iPSL.length);		
		[ribbonCoordsUnsmoothed, ~, ~, ~] = ExpandPSLs2RibbonsSim(iPSL, lineWidthTube, [6 7 8], iColor, 0);		
		ribbonCoordsUnsmoothed = reshape(ribbonCoordsUnsmoothed', numel(ribbonCoordsUnsmoothed), 1)';
		[ribbonCoordsSmoothed, ~, ~, ~] = ExpandPSLs2RibbonsSim(iPSL, lineWidthTube, [6 7 8], iColor, 1);
		ribbonCoordsSmoothed = reshape(ribbonCoordsSmoothed', numel(ribbonCoordsSmoothed), 1)';
		stressScalarFields4ColorCoding = [iPSL.principalStressList(:,1) iPSL.vonMisesStressList iPSL.cartesianStressList];
		
		fprintf(fid, '%d %.6f %.6f %.6f %.6f', [iPSL.length minorHierarchy_(ii,:)]);		
		fprintf(fid, ' %d %.6f %.6f %.6f\n', [minorPSLsGlobalAppearanceOrder(ii) iSeed]);					
		fprintf(fid, outPutFormat, pslCoords); fprintf(fid, '\n');
		fprintf(fid, outPutFormat, ribbonCoordsUnsmoothed); fprintf(fid, '\n');
		fprintf(fid, outPutFormat, ribbonCoordsSmoothed); fprintf(fid, '\n');
		for jj=1:8
			fprintf(fid, outPutFormat, stressScalarFields4ColorCoding(:,jj)'); fprintf(fid, '\n');
		end
	end
	
	%%2. write outline in quad mesh
	fprintf(fid, '%s', '#Outline');
	numVtx = size(surfaceMeshNodeCoords_,1);
	if strcmp(meshType_, 'CARTESIAN_GRID')
		fprintf(fid, '  %s', 'Cartesian'); 
	else
		fprintf(fid, '  %s', 'Unstructured'); 
	end
	if strcmp(meshType_, 'CARTESIAN_GRID') || strcmp(meshType_, 'Hex')
	
	elseif strcmp(meshType_, 'Tet')
	
	else
		error('Un-supported Mesh type!');
	end
	switch meshType_
		case 'CARTESIAN_GRID'
			tmpFaces = surfaceMeshElements_(:, [1 2 3 3 4 1]);
			tmpFaces = reshape(tmpFaces', 3, 2*size(surfaceMeshElements_,1))';
		case 'Hex'
			tmpFaces = surfaceMeshElements_(:, [1 2 3 3 4 1]);
			tmpFaces = reshape(tmpFaces', 3, 2*size(surfaceMeshElements_,1))';
		case 'Tet'
			tmpFaces = surfaceMeshElements_;
	end
	numFace = size(tmpFaces,1);
	fprintf(fid, '\n');
	fprintf(fid, '%s', '#Vertices'); fprintf(fid, ' %d\n', numVtx);
	fprintf(fid, strcat(outPutFormat, outPutFormat, outPutFormat, '\n'), surfaceMeshNodeCoords_');
	fprintf(fid, '%s', '#Faces'); fprintf(fid, ' %d\n', numFace);
	fprintf(fid, '%d %d %d\n', tmpFaces'-1);
	fclose(fid);
end