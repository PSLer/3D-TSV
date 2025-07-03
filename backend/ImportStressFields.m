function ImportStressFields(fileName)
	global boundingBox_;
	global nelx_;
	global nely_;
	global nelz_;
	global carNodMapForward_;
	global voxelizedVolume_;	
	global numNodes_;
	global nodeCoords_;
	global numEles_;
	global eNodMat_;
	global eleState_;
	global nodState_;	
	global cartesianStressField_;
	global loadingCond_; 
	global fixingCond_;	
	global eleCentroidList_;
	global silhouetteStruct_;
	global meshType_;
	global eleSize_;
	global eleSizeList_;
	global surfaceMeshNodeCoords_;
	global surfaceMeshElements_;	
	global nodStruct_; 
	global eleStruct_; 
	global boundaryElements_; 
	global numNodesPerEle_;
	global numFacesPerEle_;
	global numNodesPerFace_;
	global eleFaceIndices_;
	global nodeWiseStressField_;
	%%Read mesh and cartesian stress field
	[~,~,dataType] = fileparts(fileName);
	switch dataType
		case '.carti'
			meshType_ = 'CARTESIAN_GRID';	
			fid = fopen(fileName, 'r');
			%%Mesh
			idx = 1;
			while idx
				idx = idx + 1;
				tmp = fscanf(fid, '%s', 1);
				if strcmp(tmp, 'Type:')
					type = fscanf(fid, '%s', 1);
					switch type
						case 'NODE', nodeWiseStressField_ = 1;
						case 'ELEMENT', nodeWiseStressField_ = 0; %%Element-wise Stress Data
						otherwise, error('Un-supported Stress Data!');
					end						
				end
				if strcmp(tmp, 'Resolution:'), idx=0; break; end
				if idx>100, error('Wrong Input!'); end
			end
			tmp = fscanf(fid, '%d %d %d', [1 3]);
			nelx_ = tmp(1); nely_ = tmp(2); nelz_ = tmp(3);
			tmp = fscanf(fid, '%s', 1);
			boundingBox_ = fscanf(fid, '%f %f %f', [1 3]);
			tmp = fscanf(fid, '%s', 1);
			boundingBox_(2,:) = fscanf(fid, '%f %f %f', [1 3]);		
			tmp = fscanf(fid, '%s', 1); 
			numValidEles = fscanf(fid, '%d', 1);
			tmp = fscanf(fid, '%s', 1);
			validElements = fscanf(fid, '%d', [1, numValidEles])';
            % validElements = validElements + 1;
			% if 0==min(validElements), validElements = validElements + 1; end	
			cellVolume = zeros(nelx_*nely_*nelz_,1);
			cellVolume(validElements) = 1;
			voxelizedVolume_ = reshape(cellVolume, nely_, nelx_, nelz_);	
				
			%%Stress Field
			tmp = fscanf(fid, '%s %s %s %s %d', 5);
			tmp = fscanf(fid, '%s %s', 2); numLoadedNodes = fscanf(fid, '%d', 1);
			if numLoadedNodes>0
				loadingCond_ = fscanf(fid, '%d %f %f %f', [4, numLoadedNodes])'; 
                % loadingCond_(:,1) = loadingCond_(:,1)+1;
				% if 0==min(loadingCond_(:,1)), loadingCond_(:,1) = loadingCond_(:,1)+1; end	
			else
				loadingCond_ = [];
			end
			tmp = fscanf(fid, '%s %s', 2); numFixedNodes = fscanf(fid, '%d', 1);
			if numFixedNodes>0
				fixingCond_ = fscanf(fid, '%d', [1, numFixedNodes])';
             	% fixingCond_ = fixingCond_+1;
				% if 0==min(fixingCond_), fixingCond_ = fixingCond_+1; end
			else
				fixingCond_ = [];
			end
			tmp = fscanf(fid, '%s %s', 2); numStressComponents = fscanf(fid, '%d', 1);
			cartesianStressField_ = fscanf(fid, '%e %e %e %e %e %e', [6, numStressComponents])';	
			fclose(fid);				
			%%Recover Cartesian Mesh
			RecoverCartesianMesh();
			numNodesPerEle_ = 8;
		case '.stress'
			fid = fopen(fileName, 'r');
			fgetl(fid);
			tmp = fscanf(fid, '%s %s %s', 3);
			type = fscanf(fid, '%s', 1);
			switch type
				case 'NODE', nodeWiseStressField_ = 1;
				case 'ELEMENT', nodeWiseStressField_ = 0; %%Element-wise Stress Data
				otherwise, error('Un-supported Stress Data!');
			end			
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
			numStressComponents = fscanf(fid, '%d', 1);
			cartesianStressField_ = fscanf(fid, '%e %e %e %e %e %e', [6, numStressComponents])';		
			fclose(fid);
			boundingBox_ = [min(nodeCoords_, [], 1); max(nodeCoords_, [], 1)];
			eleSize_ = max(boundingBox_(2,:)-boundingBox_(1,:))/100;		
		otherwise
			error('Unsupported Data Format!');
	end
	%%Extract and Re-organize Silhouette into Surf-mesh for Exporting
	if strcmp(meshType_, 'CARTESIAN_GRID') || strcmp(meshType_, 'Hex')
		numNodesPerEle_ = 8;
		numFacesPerEle_ = 6;
		numNodesPerFace_ = 4;
		eleFaceIndices_ = [4 3 2 1; 5 6 7 8; 1 2 6 5; 8 7 3 4; 5 8 4 1; 2 3 7 6];
		pp = [1 2 3 4];	
	else
		numNodesPerEle_ = 4;
		numFacesPerEle_ = 4;
		numNodesPerFace_ = 3;
		eleFaceIndices_ = [3 2 1;  1 2 4;  2 3 4;  3 1 4];
		pp = [1 2 2 3];	
	end
	[surfaceMeshElements_, surfaceMeshNodeCoords_, nodState_, boundaryNodes] = ExtractBoundaryInfoFromSolidMesh();
	
	%%Extract Silhouette for Visualization
	if strcmp(meshType_, 'CARTESIAN_GRID')
		[nodPosX, nodPosY, nodPosZ] = NodalizeDesignDomain([nelx_ nely_ nelz_], boundingBox_, 'inGrid');
		iSurface = isosurface(nodPosX, nodPosY, nodPosZ, reshape(carNodMapForward_, nely_+1, nelx_+1, nelz_+1), 0);
		iCap = isocaps(nodPosX, nodPosY, nodPosZ, reshape(carNodMapForward_, nely_+1, nelx_+1, nelz_+1), 0);
		iCap.faces = size(iSurface.vertices,1) + iCap.faces;
		silhouetteStruct_.vertices = [iSurface.vertices; iCap.vertices];
		silhouetteStruct_.faces = [iSurface.faces; iCap.faces];
		silhouetteStruct_ = CompactPatchVertices(silhouetteStruct_);
	else
		silhouetteStruct_.vertices = surfaceMeshNodeCoords_;
		silhouetteStruct_.faces = surfaceMeshElements_;
	end
	
	%%element centroids
	eleNodCoordListX = nodeCoords_(:,1); eleNodCoordListX = eleNodCoordListX(eNodMat_);
	eleNodCoordListY = nodeCoords_(:,2); eleNodCoordListY = eleNodCoordListY(eNodMat_);
	eleNodCoordListZ = nodeCoords_(:,3); eleNodCoordListZ = eleNodCoordListZ(eNodMat_);
	eleCentroidList_ = [sum(eleNodCoordListX,2) sum(eleNodCoordListY,2) sum(eleNodCoordListZ,2)]/numNodesPerEle_;	
	
	%%Build Element Tree for Unstructured Hex-Mesh
	if ~strcmp(meshType_, 'CARTESIAN_GRID')
		iNodStruct = struct('adjacentEles', []); 
		nodStruct_ = repmat(iNodStruct, numNodes_, 1);
		for ii=1:numEles_
			for jj=1:numNodesPerEle_
				nodStruct_(eNodMat_(ii,jj)).adjacentEles(1,end+1) = ii;
			end
		end		
		boundaryElements_ = unique([nodStruct_(boundaryNodes).adjacentEles]);
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
end

function RecoverCartesianMesh()	
	global nelx_; 
	global nely_; 
	global nelz_; 
	global voxelizedVolume_;
	global boundingBox_;
	global numEles_; 
	global numNodes_; 
	global eleSize_;
	global nodeCoords_; 
	global eNodMat_; 
	global carEleMapBack_; 
	global carEleMapForward_;
	global carNodMapBack_; 
	global carNodMapForward_;	
	global loadingCond_; 
	global fixingCond_;
	global nodState_;
	global boundaryElements_;
	global eleCentroidList_;
	%    z
	%    |__ x
	%   / 
	%  -y                            
	%            8--------------7      	
	%			/ |			   /|	
	%          5-------------6	|
	%          |  |          |  |
	%          |  |          |  |	
	%          |  |          |  |   
	%          |  4----------|--3  
	%     	   | /           | /
	%          1-------------2             
	%			Hexahedral element
	eleSize_ = min((boundingBox_(2,:)-boundingBox_(1,:))./[nelx_ nely_ nelz_]);
	carEleMapBack_ = find(1==voxelizedVolume_);
	carEleMapBack_ = int32(carEleMapBack_);			
	numEles_ = length(carEleMapBack_);		
	carEleMapForward_ = zeros(nelx_*nely_*nelz_,1,'int32');	
	carEleMapForward_(carEleMapBack_) = (1:numEles_)';	
	nodenrs = reshape(1:(nelx_+1)*(nely_+1)*(nelz_+1), 1+nely_, 1+nelx_, 1+nelz_); nodenrs = int32(nodenrs);
	eNodVec = reshape(nodenrs(1:end-1,1:end-1,1:end-1)+1, nelx_*nely_*nelz_, 1);
	eNodMat_ = repmat(eNodVec(carEleMapBack_),1,8);
	tmp = [0 nely_+[1 0] -1 (nely_+1)*(nelx_+1)+[0 nely_+[1 0] -1]]; tmp = int32(tmp);
	for ii=1:8
		eNodMat_(:,ii) = eNodMat_(:,ii) + repmat(tmp(ii), numEles_,1);
	end
	carNodMapBack_ = unique(eNodMat_);
	numNodes_ = length(carNodMapBack_);
	carNodMapForward_ = zeros((nelx_+1)*(nely_+1)*(nelz_+1),1,'int32');
	carNodMapForward_(carNodMapBack_) = (1:numNodes_)';		
	for ii=1:8
		eNodMat_(:,ii) = carNodMapForward_(eNodMat_(:,ii));
	end
	nodeCoords_ = zeros((nelx_+1)*(nely_+1)*(nelz_+1),3);
	[nodeCoords_(:,1), nodeCoords_(:,2), nodeCoords_(:,3)] = NodalizeDesignDomain([nelx_ nely_ nelz_], boundingBox_);		
	nodeCoords_ = nodeCoords_(carNodMapBack_,:);
	if ~isempty(loadingCond_)
		loadingCond_(:,1) = carNodMapForward_(loadingCond_(:,1));
	end
	if ~isempty(fixingCond_)
		fixingCond_ = double(carNodMapForward_(fixingCond_));
	end

	numNod2ElesVec = zeros(numNodes_,1);
	for ii=1:8
		iNodes = eNodMat_(:,ii);
		numNod2ElesVec(iNodes) = numNod2ElesVec(iNodes) + 1;
	end
	nodesOutline = find(numNod2ElesVec<8);	
	nodState_ = zeros(numNodes_,1); 
	nodState_(nodesOutline) = 1;
	
	allNodes = zeros(numNodes_,1,'int32');
	allNodes(nodesOutline) = 1;	
	tmp = zeros(numEles_,1,'int32');
	for ii=1:8
		tmp = tmp + allNodes(eNodMat_(:,ii));
	end
	boundaryElements_ = int32(find(tmp>0));			
		
	%% element centroids
	blockIndex = Solving_MissionPartition(numEles_, 5.0e6);
	eleCentroidList_ = zeros(numEles_,3);
	tmp = nodeCoords_(:,1);
	for ii=1:size(blockIndex,1)	
		iSelEleNodes = eNodMat_(blockIndex(ii,1):blockIndex(ii,2),:);
		eleCentX = tmp(iSelEleNodes);
		eleCentroidList_(blockIndex(ii,1):blockIndex(ii,2),1) = sum(eleCentX,2)/8;
	end
	tmp = nodeCoords_(:,2);
	for ii=1:size(blockIndex,1)	
		iSelEleNodes = eNodMat_(blockIndex(ii,1):blockIndex(ii,2),:);
		eleCentY = tmp(iSelEleNodes);
		eleCentroidList_(blockIndex(ii,1):blockIndex(ii,2),2) = sum(eleCentY,2)/8;
	end
	tmp = nodeCoords_(:,3);
	for ii=1:size(blockIndex,1)	
		iSelEleNodes = eNodMat_(blockIndex(ii,1):blockIndex(ii,2),:);
		eleCentZ = tmp(iSelEleNodes);
		eleCentroidList_(blockIndex(ii,1):blockIndex(ii,2),3) = sum(eleCentZ,2)/8;
	end
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

function oPatchs = CompactPatchVertices(iPatchs)
	oPatchs = iPatchs;
	numOriginalVertices = size(iPatchs.vertices,1);
	numOriginalFaces = size(iPatchs.faces,1);
	validVertices = unique(iPatchs.faces);
	numValidVertices = size(validVertices,1);
	if numOriginalFaces==numOriginalVertices, return; end
	mapVerticesValid2Original = zeros(numOriginalVertices,1);
	mapVerticesValid2Original(validVertices) = (1:numValidVertices)';
	oPatchs.vertices = oPatchs.vertices(validVertices,:);
	oPatchs.faces = mapVerticesValid2Original(iPatchs.faces);
end

function blockIndex = Solving_MissionPartition(totalSize, blockSize)
	numBlocks = ceil(totalSize/blockSize);		
	blockIndex = ones(numBlocks,2);
	blockIndex(1:numBlocks-1,2) = (1:1:numBlocks-1)' * blockSize;
	blockIndex(2:numBlocks,1) = blockIndex(2:numBlocks,1) + blockIndex(1:numBlocks-1,2);
	blockIndex(numBlocks,2) = totalSize;	
end