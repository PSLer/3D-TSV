function hd = DrawMesh(varargin)
	% global silhouetteStruct_;
	global surfaceMeshNodeCoords_;
	global surfaceMeshElements_;
	
	if 0==nargin, axHandle_ = gca; else, axHandle_ = varargin{1}; end

	FV.vertices = surfaceMeshNodeCoords_;
	FV.faces = surfaceMeshElements_;
	hd = patch(axHandle_, FV); 
	hold(axHandle_, 'on');
	set(hd, 'FaceColor', [65 174 118]/255, 'FaceAlpha', 1.0, 'EdgeColor', 'k');
	
	view(axHandle_, 3);
	camproj(axHandle_, 'perspective');
	axis(axHandle_, 'equal'); 
	axis(axHandle_, 'tight');
	axis(axHandle_, 'off');
	
	%%Lighting, Reflection
	lighting(axHandle_, 'gouraud');
	material(axHandle_, 'dull');
	camlight(axHandle_, 'headlight', 'infinite');
end