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
	global surfaceQuadMeshNodeCoords_;
	global surfaceQuadMeshElements_;	
	global nodStruct_; 
	global eleStruct_; 
	global boundaryElements_; 
	
	%%Read mesh and cartesian stress field
	[~,~,dataType] = fileparts(fileName);
	switch dataType
		case '.carti'
			meshType_ = 'CARTESIAN_GRID';
			[nelx_, nely_, nelz_, boundingBox_, voxelizedVolume_, loadingCond_, fixingCond_, cartesianStressField_] = ...
				ReadDataSimulatedOnCartesianMesh_carti(fileName);
			%%Recover Cartesian Mesh
			RecoverCartesianMesh();
		case '.stress'
			meshType_ = 'UNSTRUCTURED_GRID';
			[numNodes_, nodeCoords_, numEles_, eNodMat_, loadingCond_, fixingCond_, cartesianStressField_] = ...
				ReadDataSimulatedOnUnstructuredHexMesh_stress(fileName);
			%%Extract Boundary Element Info.
			faceIndex = eNodMat_(:, [4 3 2 1  5 6 7 8  1 2 6 5  8 7 3 4  5 8 4 1  2 3 7 6])';
			faceIndex = reshape(faceIndex(:), 4, 6*numEles_);
			tmp = sort(faceIndex,1)';
			[~, ia, ic] = unique(tmp, 'rows');
			numRawPatchs = 6*numEles_;
			patchState = zeros(length(ia),1);
			for ii=1:numRawPatchs
				patchState(ic(ii)) = patchState(ic(ii)) + 1;
			end
			patchIndexOnBoundary = ia(1==patchState);
			boundaryPatchs = faceIndex(:,patchIndexOnBoundary');
			nodesOutline = unique(boundaryPatchs);
			nodState_ = zeros(numNodes_,1);
			nodState_(nodesOutline) = 1;					
		otherwise
			error('Unsupported Data Format!');
	end
	
	%%Extract and Re-organize Silhouette into Quad-mesh for Exporting
	faceIndex = eNodMat_(:, [4 3 2 1  5 6 7 8  1 2 6 5  8 7 3 4  5 8 4 1  2 3 7 6])';
	faceIndex = reshape(faceIndex(:), 4, 6*numEles_);	
	tmp = nodState_(faceIndex'); 
	tmp = sum(tmp,2);
	boundaryEleFace = faceIndex(:,find(4==tmp)');		
	boundaryNode = find(1==nodState_);
	surfaceQuadMeshNodeCoords_ = nodeCoords_(boundaryNode,:);
	tmp = zeros(numNodes_,1); 
	tmp(boundaryNode) = (1:length(boundaryNode))';
	surfaceQuadMeshElements_ = tmp(boundaryEleFace');
	
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
		silhouetteStruct_.vertices = surfaceQuadMeshNodeCoords_;
		silhouetteStruct_.faces = surfaceQuadMeshElements_;
	end
	
	%%element centroids
	eleNodCoordListX = nodeCoords_(:,1); eleNodCoordListX = eleNodCoordListX(eNodMat_);
	eleNodCoordListY = nodeCoords_(:,2); eleNodCoordListY = eleNodCoordListY(eNodMat_);
	eleNodCoordListZ = nodeCoords_(:,3); eleNodCoordListZ = eleNodCoordListZ(eNodMat_);
	eleCentroidList_ = [sum(eleNodCoordListX,2) sum(eleNodCoordListY,2) sum(eleNodCoordListZ,2)]/8;	
	
	%%Build Element Tree for Unstructured Hex-Mesh
	if strcmp(meshType_, 'UNSTRUCTURED_GRID')
		boundingBox_ = [min(nodeCoords_, [], 1); max(nodeCoords_, [], 1)];
		eleSize_ = max(boundingBox_(2,:)-boundingBox_(1,:))/100;
		%% build element three		
		iNodStruct = struct('adjacentEles', []); 
		nodStruct_ = repmat(iNodStruct, numNodes_, 1);
		for ii=1:numEles_
			for jj=1:8
				nodStruct_(eNodMat_(ii,jj)).adjacentEles(1,end+1) = ii;
			end
		end		
		boundaryElements_ = unique([nodStruct_(boundaryNode).adjacentEles]);
		boundaryElements_ = boundaryElements_(:);
		eleState_ = zeros(numEles_,1);
		eleState_(boundaryElements_,1) = 1;
		eleFaces = [4 3 2 1; 5 6 7 8; 1 2 6 5; 8 7 3 4; 5 8 4 1; 2 3 7 6];
		iEleStruct = struct('faceCentres', [], 'faceNormals', []); %%pure-Hex
		eleStruct_ = repmat(iEleStruct, numEles_, 1);
		for ii=1:numEles_
			iNodes = eNodMat_(ii,:);
			iEleVertices = nodeCoords_(iNodes, :);
			iEleFacesX = iEleVertices(:,1); iEleFacesX = iEleFacesX(eleFaces);
			iEleFacesY = iEleVertices(:,2); iEleFacesY = iEleFacesY(eleFaces);
			iEleFacesZ = iEleVertices(:,3); iEleFacesZ = iEleFacesZ(eleFaces);				
			ACs = [iEleFacesX(:,1)-iEleFacesX(:,3) iEleFacesY(:,1)-iEleFacesY(:,3) iEleFacesZ(:,1)-iEleFacesZ(:,3)];
			BDs = [iEleFacesX(:,2)-iEleFacesX(:,4) iEleFacesY(:,2)-iEleFacesY(:,4) iEleFacesZ(:,2)-iEleFacesZ(:,4)];
			iACxBD = cross(ACs,BDs); 
			aveNormal = iACxBD ./ vecnorm(iACxBD,2,2);			
			tmp = iEleStruct;			
			%% tmp.faceNormals = aveNormal;
			%% in case the node orderings on each element face are not constant
			tmp.faceCentres = [sum(iEleFacesX,2) sum(iEleFacesY,2) sum(iEleFacesZ,2)]/4;
			iEleCt = eleCentroidList_(ii,:);
			refVecs = iEleCt - tmp.faceCentres; refVecs = refVecs ./ vecnorm(refVecs,2,2);
			dirEval = acos(sum(refVecs .* aveNormal, 2));
			dirDes = ones(6,1); dirDes(dirEval<pi/2) = -1;
			faceNormals = dirDes .* aveNormal;
			tmp.faceNormals = faceNormals;
			eleStruct_(ii) = tmp;
		end
		
		%% Evaluate Element Sizes
		tmpSizeList = zeros(6,numEles_);
		for ii=1:numEles_
			tmpSizeList(:,ii) = vecnorm(eleCentroidList_(ii,:)-eleStruct_(ii).faceCentres,2,2);
		end
		eleSizeList_ = 2*min(tmpSizeList,[],1)';		
	end
end

function [nx, ny, nz, boundingBox, cellVolume, loadingCond, fixingCond, cartesianStressField] = ...
				ReadDataSimulatedOnCartesianMesh_carti(fileName)
	fid = fopen(fileName, 'r');
	%%Mesh
	idx = 1;
	while idx
		idx = idx + 1;
		tmp = fscanf(fid, '%s', 1);
		if strcmp(tmp, 'Resolution:'), idx=0; break; end
		if idx>100, error('Wrong Input!'); end
	end
	tmp = fscanf(fid, '%d %d %d', [1 3]);
	nx = tmp(1); ny = tmp(2); nz = tmp(3);
	tmp = fscanf(fid, '%s', 1);
	boundingBox = fscanf(fid, '%f %f %f', [1 3]);
	tmp = fscanf(fid, '%s', 1);
	boundingBox(2,:) = fscanf(fid, '%f %f %f', [1 3]);		
	tmp = fscanf(fid, '%s', 1); 
	numValidEles = fscanf(fid, '%d', 1);
	tmp = fscanf(fid, '%s', 1);
	validElements = fscanf(fid, '%d', [1, numValidEles])';
	validElements = validElements + 1;
	cellVolume = zeros(nx*ny*nz,1);
	cellVolume(validElements) = 1;
	cellVolume = reshape(cellVolume, ny, nx, nz);	
		
	%%Stress Field
	tmp = fscanf(fid, '%s %s %s %s %d', 5);
	tmp = fscanf(fid, '%s %s', 2); numLoadedNodes = fscanf(fid, '%d', 1);
	if numLoadedNodes>0
		tmp = fscanf(fid, '%d %f %f %f', [4, numLoadedNodes]); 
		tmp(1,:) = tmp(1,:)+1; 
		loadingCond = tmp';
	else
		loadingCond = [];
	end
	tmp = fscanf(fid, '%s %s', 2); numFixedNodes = fscanf(fid, '%d', 1);
	if numFixedNodes>0
		tmp = fscanf(fid, '%d', [1, numFixedNodes]); 
		fixingCond = tmp'+1;
	else
		fixingCond = [];
	end
	tmp = fscanf(fid, '%s %s', 2); numValidNods = fscanf(fid, '%d', 1);
	cartesianStressField = fscanf(fid, '%e %e %e %e %e %e', [6, numValidNods])';	
	fclose(fid);
end

function [numNodes, nodeCoords, numEles, eNodMat, loadingCond, fixingCond, cartesianStressField] = ...
				ReadDataSimulatedOnUnstructuredHexMesh_stress(fileName)
	fid = fopen(fileName, 'r');
	%%Mesh
	idx = 1;
	while idx
		idx = idx + 1;
		tmp = fscanf(fid, '%s', 1);
		if strcmp(tmp, 'Vertices:'), idx=0; break; end
		if idx>100, error('Wrong Input!'); end
	end
	numNodes = fscanf(fid, '%d', 1);
	nodeCoords = fscanf(fid, '%f %f %f', [3, numNodes])'; 
	
	tmp = fscanf(fid, '%s', 1);
	numEles = fscanf(fid, '%d', 1);
	eNodMat = fscanf(fid, '%d %d %d %d %d %d %d %d', [8, numEles])'; 
	
	%%Stress Field
	tmp = fscanf(fid, '%s %s', 2); 
	numLoadedNodes = fscanf(fid, '%d', 1);
	if numLoadedNodes>0
		tmp = fscanf(fid, '%d %e %e %e', [4, numLoadedNodes]); 
		loadingCond = tmp';
	else
		loadingCond = [];
	end
	tmp = fscanf(fid, '%s %s', 2); 
	numFixedNodes = fscanf(fid, '%d', 1);
	if numFixedNodes>0
		tmp = fscanf(fid, '%d', [1, numFixedNodes]); 
		fixingCond = tmp(:);
	else
		fixingCond = [];
	end
	tmp = fscanf(fid, '%s %s', 2); 
	numValidNods = fscanf(fid, '%d', 1);
	cartesianStressField = fscanf(fid, '%e %e %e %e %e %e', [6, numValidNods])';	
	fclose(fid);
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
	eleNodCoordListX = nodeCoords_(:,1); eleNodCoordListX = eleNodCoordListX(eNodMat_);
	eleNodCoordListY = nodeCoords_(:,2); eleNodCoordListY = eleNodCoordListY(eNodMat_);
	eleNodCoordListZ = nodeCoords_(:,3); eleNodCoordListZ = eleNodCoordListZ(eNodMat_);
	eleCentroidList_ = [sum(eleNodCoordListX,2) sum(eleNodCoordListY,2) sum(eleNodCoordListZ,2)]/8;	
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