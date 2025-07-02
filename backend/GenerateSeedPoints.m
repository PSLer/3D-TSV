function GenerateSeedPoints(seedStrategy, seedDensCtrl)
	global meshType_;
	global nodeCoords_;
	global nodState_;
	global loadingCond_;
	global fixingCond_;
	global seedPointsHistory_;
	global seedAssociatedEles_;
	global carNodMapBack_;
	global nelx_; 
	global nely_; 
	global nelz_;
	global nodStruct_;
	global boundaryElements_;
	global eleCentroidList_;
	
	if 0>=seedDensCtrl, seedDensCtrl = 4; end %% in case 0>=minimumEpsilon and volumeSeedingOpt is not defined
	switch seedStrategy
		case 'Volume'
			if strcmp(meshType_, 'CARTESIAN_GRID')
				if seedDensCtrl > min([nelx_ nely_ nelz_])/3, seedDensCtrl = 4; end %% in case seedDensCtrl is wrongly selected
				validNodesVolume = zeros((nelx_+1)*(nely_+1)*(nelz_+1),1);
				validNodesVolume(carNodMapBack_) = (1:length(carNodMapBack_))';
				validNodesVolume = reshape(validNodesVolume, nely_+1, nelx_+1, nelz_+1);
				sampledNodes = validNodesVolume(seedDensCtrl+1:seedDensCtrl:nely_+1-seedDensCtrl, ...
					seedDensCtrl+1:seedDensCtrl:nelx_+1-seedDensCtrl, seedDensCtrl+1:seedDensCtrl:nelz_+1-seedDensCtrl);	
				sampledNodes = reshape(sampledNodes, numel(sampledNodes), 1);
				sampledNodes(0==sampledNodes) = [];
				nodesOnBoundary = find(1==nodState_);
				sampledNodes = setdiff(sampledNodes, nodesOnBoundary);
				seedPointsHistory_ = nodeCoords_(sampledNodes,:);
				seedAssociatedEles_ = ones(size(seedPointsHistory_,1),1); %%unused
			else				
				tarNodes = find(0==nodState_);
				tarNodes = tarNodes(1:seedDensCtrl:end,1);
				seedPointsHistory_ = nodeCoords_(tarNodes,:); %%exclude the boundary nodes
				seedAssociatedEles_ = zeros(size(tarNodes));
				for ii=1:length(tarNodes)
					iNode = tarNodes(ii);
					seedAssociatedEles_(ii) = nodStruct_(iNode).adjacentEles(1);
				end
			end
		case 'Surface'
			potentialElements = boundaryElements_(1:seedDensCtrl:end,1);
			seedPointsHistory_ = eleCentroidList_(potentialElements,:);
			seedAssociatedEles_ = potentialElements; %%unused for Cartesian mesh
		case 'LoadingArea'			
			if isempty(loadingCond_)
				seedPointsHistory_ = [];
				warning('There is no Loaded Node Available!'); return;
			end
			if strcmp(meshType_, 'CARTESIAN_GRID')
				seedPointsHistory_ = nodeCoords_(loadingCond_(1:seedDensCtrl:end,1),:);
				seedAssociatedEles_ = ones(size(seedPointsHistory_,1),1); %%unused
			else
				potentialElements = unique([nodStruct_(loadingCond_(:,1)).adjacentEles]);
				potentialElements = potentialElements(:);
				potentialElements = potentialElements(1:seedDensCtrl:end,1);
				seedPointsHistory_ = eleCentroidList_(potentialElements,:);
				seedAssociatedEles_ = potentialElements;
			end	
		case 'FixedArea' 
			if isempty(fixingCond_)
				seedPointsHistory_ = [];
				warning('There is no Fixed Node Available!'); return;		
			end
			if strcmp(meshType_, 'CARTESIAN_GRID')
				seedPointsHistory_ = nodeCoords_(fixingCond_(1:seedDensCtrl:end,1),:);
				seedAssociatedEles_ = ones(size(seedPointsHistory_,1),1); %%unused
			else				
				potentialElements = unique([nodStruct_(fixingCond_).adjacentEles]);
				potentialElements = potentialElements(:);
				potentialElements = potentialElements(1:seedDensCtrl:end,1);
				seedPointsHistory_ = eleCentroidList_(potentialElements,:);
				seedAssociatedEles_ = potentialElements;
			end
		otherwise
			error('Unsupported Seed Type!');
	end
end