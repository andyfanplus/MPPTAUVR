#include"vtkObjectFactory.h"
#include"vtkPRenderingVolumeOpenGL2ObjectFactory.h"
#include"vtkOpenGLPProjectedTetrahedraMapper.h"

VTK_CREATE_CREATE_FUNCTION(vtkOpenGLPProjectedTetrahedraMapper);


vtkPRenderingVolumeOpenGL2ObjectFactory::vtkPRenderingVolumeOpenGL2ObjectFactory()
{
	this->RegisterOverride("vtkPProjectedTetrahedraMapper",
		"vtkOpenGLPProjectedTetrahedraMapper",
		"test vertex factory override",
		1,
		vtkObjectFactoryCreatevtkOpenGLPProjectedTetrahedraMapper);
}