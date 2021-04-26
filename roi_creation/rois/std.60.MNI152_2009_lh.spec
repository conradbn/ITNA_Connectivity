# MapIcosahedron generated spec file
#History: [taylorpa3@cn0686: Fri May 25 00:40:45 2018] MapIcosahedron -spec MNI152_2009_lh.spec -ld 60 -dset_map lh.thickness.gii.dset -dset_map lh.curv.gii.dset -dset_map lh.sulc.gii.dset -prefix std.60.

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
	SurfaceName = ./std.60.lh.smoothwm.gii
	LocalDomainParent = ./SAME
	LabelDset = ./std.60.lh.aparc.a2009s.annot.niml.dset
	SurfaceState = std.smoothwm
	EmbedDimension = 3
	Anatomical = Y
	LocalCurvatureParent = ./SAME

NewSurface
	SurfaceFormat = ASCII
	SurfaceType = GIFTI
	SurfaceName = ./std.60.lh.pial.gii
	LocalDomainParent = ./std.60.lh.smoothwm.gii
	SurfaceState = std.pial
	EmbedDimension = 3
	Anatomical = Y
	LocalCurvatureParent = ./std.60.lh.smoothwm.gii

NewSurface
	SurfaceFormat = ASCII
	SurfaceType = GIFTI
	SurfaceName = ./std.60.lh.inflated.gii
	LocalDomainParent = ./std.60.lh.smoothwm.gii
	SurfaceState = std.inflated
	EmbedDimension = 3
	Anatomical = N
	LocalCurvatureParent = ./std.60.lh.smoothwm.gii

NewSurface
	SurfaceFormat = ASCII
	SurfaceType = GIFTI
	SurfaceName = ./std.60.lh.sphere.gii
	LocalDomainParent = ./std.60.lh.smoothwm.gii
	SurfaceState = std.sphere
	EmbedDimension = 3
	Anatomical = N
	LocalCurvatureParent = ./std.60.lh.smoothwm.gii

NewSurface
	SurfaceFormat = ASCII
	SurfaceType = GIFTI
	SurfaceName = ./std.60.lh.white.gii
	LocalDomainParent = ./std.60.lh.smoothwm.gii
	SurfaceState = std.white
	EmbedDimension = 3
	Anatomical = Y
	LocalCurvatureParent = ./std.60.lh.smoothwm.gii

NewSurface
	SurfaceFormat = ASCII
	SurfaceType = GIFTI
	SurfaceName = ./std.60.lh.sphere.reg.gii
	LocalDomainParent = ./std.60.lh.smoothwm.gii
	SurfaceState = std.sphere.reg
	EmbedDimension = 3
	Anatomical = N
	LocalCurvatureParent = ./std.60.lh.smoothwm.gii

NewSurface
	SurfaceFormat = ASCII
	SurfaceType = GIFTI
	SurfaceName = ./std.60.lh.inf_200.gii
	LocalDomainParent = ./std.60.lh.smoothwm.gii
	SurfaceState = std.inf_200
	EmbedDimension = 3
	Anatomical = N
	LocalCurvatureParent = ./std.60.lh.smoothwm.gii
