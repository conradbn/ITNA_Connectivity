# MapIcosahedron generated spec file
#History: [taylorpa3@cn0686: Fri May 25 00:39:58 2018] MapIcosahedron -spec MNI152_2009_rh.spec -ld 141 -dset_map rh.thickness.gii.dset -dset_map rh.curv.gii.dset -dset_map rh.sulc.gii.dset -prefix std.141.

#define the group
	Group = MNI152_2009

#define various States
	StateDef = std.smoothwm
	StateDef = std.pial
	StateDef = std.inflated
	StateDef = std.sphere
	StateDef = std.white
	StateDef = std.sphere.reg
	StateDef = std.inf_200

NewSurface
	SurfaceFormat = ASCII
	SurfaceType = GIFTI
	SurfaceName = ./std.141.rh.smoothwm.gii
	LocalDomainParent = ./SAME
	LabelDset = ./std.141.rh.aparc.a2009s.annot.niml.dset
	SurfaceState = std.smoothwm
	EmbedDimension = 3
	Anatomical = Y
	LocalCurvatureParent = ./SAME

NewSurface
	SurfaceFormat = ASCII
	SurfaceType = GIFTI
	SurfaceName = ./std.141.rh.pial.gii
	LocalDomainParent = ./std.141.rh.smoothwm.gii
	SurfaceState = std.pial
	EmbedDimension = 3
	Anatomical = Y
	LocalCurvatureParent = ./std.141.rh.smoothwm.gii

NewSurface
	SurfaceFormat = ASCII
	SurfaceType = GIFTI
	SurfaceName = ./std.141.rh.inflated.gii
	LocalDomainParent = ./std.141.rh.smoothwm.gii
	SurfaceState = std.inflated
	EmbedDimension = 3
	Anatomical = N
	LocalCurvatureParent = ./std.141.rh.smoothwm.gii

NewSurface
	SurfaceFormat = ASCII
	SurfaceType = GIFTI
	SurfaceName = ./std.141.rh.sphere.gii
	LocalDomainParent = ./std.141.rh.smoothwm.gii
	SurfaceState = std.sphere
	EmbedDimension = 3
	Anatomical = N
	LocalCurvatureParent = ./std.141.rh.smoothwm.gii

NewSurface
	SurfaceFormat = ASCII
	SurfaceType = GIFTI
	SurfaceName = ./std.141.rh.white.gii
	LocalDomainParent = ./std.141.rh.smoothwm.gii
	SurfaceState = std.white
	EmbedDimension = 3
	Anatomical = Y
	LocalCurvatureParent = ./std.141.rh.smoothwm.gii

NewSurface
	SurfaceFormat = ASCII
	SurfaceType = GIFTI
	SurfaceName = ./std.141.rh.sphere.reg.gii
	LocalDomainParent = ./std.141.rh.smoothwm.gii
	SurfaceState = std.sphere.reg
	EmbedDimension = 3
	Anatomical = N
	LocalCurvatureParent = ./std.141.rh.smoothwm.gii

NewSurface
	SurfaceFormat = ASCII
	SurfaceType = GIFTI
	SurfaceName = ./std.141.rh.inf_200.gii
	LocalDomainParent = ./std.141.rh.smoothwm.gii
	SurfaceState = std.inf_200
	EmbedDimension = 3
	Anatomical = N
	LocalCurvatureParent = ./std.141.rh.smoothwm.gii
