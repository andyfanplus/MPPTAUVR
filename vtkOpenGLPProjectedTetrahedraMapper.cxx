#include "vtkOpenGLPProjectedTetrahedraMapper.h"
#include "vtkCamera.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellIterator.h"
#include "vtkFloatArray.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkMath.h"
#include "vtkMatrix3x3.h"
#include "vtkMatrix4x4.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkOpenGLCamera.h"
#include "vtkOpenGLFramebufferObject.h"
#include "vtkOpenGLIndexBufferObject.h"
#include "vtkOpenGLRenderWindow.h"
#include "vtkOpenGLShaderCache.h"
#include "vtkOpenGLVertexArrayObject.h"
#include "vtkOpenGLVertexBufferObject.h"
#include "vtkPointData.h"
#include "vtkRenderer.h"
#include "vtkShaderProgram.h"
#include "vtkSmartPointer.h"
#include "vtkTimerLog.h"
#include "vtkUnsignedCharArray.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPVisibilitySort.h"//PARALLEL
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"

#include "vtkOpenGLError.h"

#include <cmath>
#include <algorithm>
#include<omp.h>


const char *vtkglProjectedTetrahedraVS1 =
"//VTK::System::Dec\n"
"\n"
"/*=========================================================================\n"
"\n"
"  Program:   Visualization Toolkit\n"
"  Module:    vtkglProjectedTetrahedra.glsl\n"
"\n"
"  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen\n"
"  All rights reserved.\n"
"  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.\n"
"\n"
"     This software is distributed WITHOUT ANY WARRANTY; without even\n"
"     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR\n"
"     PURPOSE.  See the above copyright notice for more information.\n"
"\n"
"=========================================================================*/\n"
"\n"
"\n"
"// all variables that represent positions or directions have a suffix\n"
"// indicating the coordinate system they are in. The possible values are\n"
"// MC - Model Coordinates\n"
"// WC - WC world coordinates\n"
"// VC - View Coordinates\n"
"// DC - Display Coordinates\n"
"in vec4 vertexDC;\n"
"in vec3 scalarColor;\n"
"in float depthArray;\n"
"in float attenuationArray;\n"
"\n"
"out float fdepth;\n"
"out float fattenuation;\n"
"out vec3 fcolor;\n"
"\n"
"void main()\n"
"{\n"
"  fcolor = scalarColor;\n"
"  fdepth = depthArray;\n"
"  fattenuation = attenuationArray;\n"
"  gl_Position = vertexDC;\n"
"}\n"
"";

const char *vtkglProjectedTetrahedraFS1 =
"//VTK::System::Dec\n"
"\n"
"/*=========================================================================\n"
"\n"
"  Program:   Visualization Toolkit\n"
"  Module:    vtkglprojectedTetrahdraFS.glsl\n"
"\n"
"  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen\n"
"  All rights reserved.\n"
"  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.\n"
"\n"
"     This software is distributed WITHOUT ANY WARRANTY; without even\n"
"     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR\n"
"     PURPOSE.  See the above copyright notice for more information.\n"
"\n"
"=========================================================================*/\n"
"\n"
"//VTK::Output::Dec\n"
"\n"
"varying vec3 fcolor;\n"
"varying float fdepth;\n"
"varying float fattenuation;\n"
"\n"
"void main()\n"
"{\n"
"  // the following exp is done in the fragment shader\n"
"  // because linear interpolation (from the VS) of the resulting\n"
"  // value would not match the exp of the interpolated\n"
"  // source values\n"
"  float opacity = 1.0 - exp(-1.0*fattenuation*fdepth);\n"
"\n"
"\n"
"  gl_FragData[0] =  vec4(fcolor,opacity);\n"
"\n"
"  if (gl_FragData[0].a <= 0.0)\n"
"    {\n"
"    discard;\n"
"    }\n"
"}\n"
"\n";

static int tet_edges[6][2] = { { 0,1 },{ 1,2 },{ 2,0 },
{ 0,3 },{ 1,3 },{ 2,3 } };

const int SqrtTableSize = 2048;


//-----------------------------------------------------------------------------
vtkStandardNewMacro(vtkOpenGLPProjectedTetrahedraMapper);


void vtkOpenGLPProjectedTetrahedraMapper::SetNumberOfVST(int vst)
{
	this->VST = vst;
}
void vtkOpenGLPProjectedTetrahedraMapper::SetNumberOfCDT(int cdt)
{
	this->CDT = cdt;
}



//-----------------------------------------------------------------------------
vtkOpenGLPProjectedTetrahedraMapper::vtkOpenGLPProjectedTetrahedraMapper()
{
	this->TransformedPoints = vtkFloatArray::New();
	this->Colors = vtkUnsignedCharArray::New();
	this->LastProperty = NULL;
	this->MaxCellSize = 0;
	this->GaveError = 0;
	this->SqrtTable = new float[SqrtTableSize];
	this->SqrtTableBias = 0.0;
	this->Initialized = false;
	this->CurrentFBOWidth = -1;
	this->CurrentFBOHeight = -1;
	this->FloatingPointFrameBufferResourcesAllocated = false;
	this->Framebuffer = vtkOpenGLFramebufferObject::New();
	this->UseFloatingPointFrameBuffer = true;
	this->CanDoFloatingPointFrameBuffer = false;
	this->HasHardwareSupport = false;
	this->VBO = vtkOpenGLVertexBufferObject::New();

	//this->packedVBO = new float[45000000];
	//this->indexArray = new unsigned int[18000000];
	//this->cellArray = new unsigned int[6000000];
	this->packedVBO = new float[120000000];
	this->indexArray = new unsigned int[64000000];
	this->cellArray = new unsigned int[64000000];
}

//-----------------------------------------------------------------------------
vtkOpenGLPProjectedTetrahedraMapper::~vtkOpenGLPProjectedTetrahedraMapper()
{
	this->ReleaseGraphicsResources(NULL);
	this->TransformedPoints->Delete();
	this->Colors->Delete();
	delete[] this->SqrtTable;
	this->VBO->Delete();
	this->Framebuffer->Delete();
	delete[] this->packedVBO;
	delete[] this->indexArray;
	delete[] this->cellArray;
}

//-----------------------------------------------------------------------------
void vtkOpenGLPProjectedTetrahedraMapper::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "PVisibilitySort: " << this->VisibilitySort << endl;
	os << indent << "UseFloatingPointFrameBuffer: "
		<< (this->UseFloatingPointFrameBuffer ? "True" : "False") << endl;
}

//-----------------------------------------------------------------------------
bool vtkOpenGLPProjectedTetrahedraMapper::IsSupported(vtkRenderWindow *rwin)
{
	vtkOpenGLRenderWindow *context = vtkOpenGLRenderWindow::SafeDownCast(rwin);
	if (!context)
	{
		vtkErrorMacro(
			<< "Support for " << rwin->GetClassName() << " not implemented");
		return false;
	}

	// use render to FBO when it's supported
	this->CanDoFloatingPointFrameBuffer = false;
	if (this->UseFloatingPointFrameBuffer)
	{
		this->CanDoFloatingPointFrameBuffer = true;
	}


	return true;
}

//-----------------------------------------------------------------------------
void vtkOpenGLPProjectedTetrahedraMapper::Initialize(vtkRenderer *renderer)
{
	if (this->Initialized)
	{
		return;
	}

	this->Initialized = true;

	vtkOpenGLRenderWindow *renwin
		= vtkOpenGLRenderWindow::SafeDownCast(renderer->GetRenderWindow());
	this->HasHardwareSupport = renwin != NULL && this->IsSupported(renwin);
	if (!this->HasHardwareSupport)
	{
		// this is an error since there's no fallback.
		vtkErrorMacro("The required extensions are not supported.");
	}
}

//-----------------------------------------------------------------------------
bool vtkOpenGLPProjectedTetrahedraMapper::AllocateFOResources(vtkRenderer *r)
{
	vtkOpenGLClearErrorMacro();

	int *size = r->GetSize();

	if (this->UseFloatingPointFrameBuffer
		&& this->CanDoFloatingPointFrameBuffer
		&& (!this->FloatingPointFrameBufferResourcesAllocated
			|| (size[0] != this->CurrentFBOWidth)
			|| (size[0] != this->CurrentFBOHeight)))
	{
		vtkOpenGLRenderWindow *rw =
			static_cast<vtkOpenGLRenderWindow *>(r->GetRenderWindow());

		if (!this->FloatingPointFrameBufferResourcesAllocated)
		{
			// determine if we have MSAA
			GLint winSampleBuffers = 0;
			glGetIntegerv(GL_SAMPLE_BUFFERS, &winSampleBuffers);
			GLint winSamples = 0;
			if (winSampleBuffers)
			{
				glGetIntegerv(GL_SAMPLES, &winSamples);
			}

			int dsize = rw->GetDepthBufferSize();
			if (dsize == 0)
			{
				dsize = 24;
			}

			vtkOpenGLFramebufferObject *fo = this->Framebuffer;
			fo->SetContext(rw);
			fo->SaveCurrentBindingsAndBuffers();

			const char *desc;

			// if we failed to get a framebuffer and we wanted
			// multisamples, then try again without multisamples
			if (!fo->PopulateFramebuffer(size[0], size[1],
				true, // use textures
				1, VTK_FLOAT, // 1 color buffer of float
				true, dsize, // yes depth buffer
				winSamples) // possibly multisampled
				&& winSamples > 0)
			{
				fo->PopulateFramebuffer(size[0], size[1],
					true, // use textures
					1, VTK_FLOAT, // 1 color buffer of float
					true, dsize, // yes depth buffer
					0); // no multisamples
			}

			this->FloatingPointFrameBufferResourcesAllocated = true;

			if (!fo->GetFrameBufferStatus(fo->GetDrawMode(), desc))
			{
				vtkWarningMacro(
					"Missing FBO support. The algorithm may produce visual artifacts.");
				this->CanDoFloatingPointFrameBuffer = false;
				fo->RestorePreviousBindingsAndBuffers();
				return false;
			}
			this->Framebuffer->UnBind();
			fo->RestorePreviousBindingsAndBuffers();
			this->CanDoFloatingPointFrameBuffer = true;
		}
		else
		{
			// need resize
			vtkOpenGLFramebufferObject *fo = this->Framebuffer;
			fo->SaveCurrentBindingsAndBuffers();
			fo->Bind();
			fo->Resize(size[0], size[1]);
			this->Framebuffer->UnBind();
			fo->RestorePreviousBindingsAndBuffers();
		}
		this->CurrentFBOWidth = size[0];
		this->CurrentFBOHeight = size[1];
	}
	return true;
}

//-----------------------------------------------------------------------------
void vtkOpenGLPProjectedTetrahedraMapper::ReleaseGraphicsResources(vtkWindow *win)
{
	this->Initialized = false;

	if (this->FloatingPointFrameBufferResourcesAllocated)
	{
		this->FloatingPointFrameBufferResourcesAllocated = false;
		this->Framebuffer->ReleaseGraphicsResources(win);
	}

	this->VBO->ReleaseGraphicsResources();
	this->Tris.ReleaseGraphicsResources(win);

	this->Superclass::ReleaseGraphicsResources(win);
}

//-----------------------------------------------------------------------------
void vtkOpenGLPProjectedTetrahedraMapper::Render(vtkRenderer *renderer,
	vtkVolume *volume)
{
	vtkOpenGLClearErrorMacro();

	// load required extensions
	this->Initialize(renderer);

	if (!this->HasHardwareSupport)
	{
		return;
	}

	// make sure our shader program is loaded and ready to go
	vtkOpenGLRenderWindow *renWin =
		vtkOpenGLRenderWindow::SafeDownCast(renderer->GetRenderWindow());

	if (renWin == NULL)
	{
		vtkErrorMacro("Invalid vtkOpenGLRenderWindow");
	}

	vtkUnstructuredGridBase *input = this->GetInput();
	vtkVolumeProperty *property = volume->GetProperty();

	// has something changed that would require us to recreate the shader?
	if (!this->Tris.Program)
	{
		// build the shader source code
		std::string VSSource = vtkglProjectedTetrahedraVS1;
		std::string FSSource = vtkglProjectedTetrahedraFS1;
		std::string GSSource;

		// compile and bind it if needed
		vtkShaderProgram *newShader =
			renWin->GetShaderCache()->ReadyShaderProgram(VSSource.c_str(),
				FSSource.c_str(),
				GSSource.c_str());

		// if the shader changed reinitialize the VAO
		if (newShader != this->Tris.Program)
		{
			this->Tris.Program = newShader;
			this->Tris.VAO->ShaderProgramChanged(); // reset the VAO as the shader has changed
		}

		this->Tris.ShaderSourceTime.Modified();
	}
	else
	{
		renWin->GetShaderCache()->ReadyShaderProgram(this->Tris.Program);
	}

	// Check to see if input changed.
	if ((this->InputAnalyzedTime < this->MTime)
		|| (this->InputAnalyzedTime < input->GetMTime()))
	{
		this->GaveError = 0;
		float max_cell_size2 = 0;

		if (input->GetNumberOfCells() == 0)
		{
			// Apparently, the input has no cells.  Just do nothing.
			return;
		}
		unsigned int *tmpCell = this->cellArray;
		vtkSmartPointer<vtkCellIterator> cellIter =
			vtkSmartPointer<vtkCellIterator>::Take(input->NewCellIterator());
		for (cellIter->InitTraversal(); !cellIter->IsDoneWithTraversal();
			cellIter->GoToNextCell())
		{
			vtkIdType npts = cellIter->GetNumberOfPoints();
			if (npts != 4)
			{
				if (!this->GaveError)
				{
					vtkErrorMacro("Encountered non-tetrahedra cell!");
					this->GaveError = 1;
				}
				continue;
			}
			vtkIdType *pts = cellIter->GetPointIds()->GetPointer(0);
			*(tmpCell++) = pts[0];
			*(tmpCell++) = pts[1];
			*(tmpCell++) = pts[2];
			*(tmpCell++) = pts[3];
			for (int j = 0; j < 6; j++)
			{
				double p1[3], p2[3];
				input->GetPoint(pts[tet_edges[j][0]], p1);
				input->GetPoint(pts[tet_edges[j][1]], p2);
				float size2 = (float)vtkMath::Distance2BetweenPoints(p1, p2);
				if (size2 > max_cell_size2)
				{
					max_cell_size2 = size2;
				}
			}
		}
		this->MaxCellSize = (float)sqrt(max_cell_size2);

		// Build a sqrt lookup table for measuring distances.  During perspective
		// modes we have to take a lot of square roots, and a table is much faster
		// than calling the sqrt function.
		this->SqrtTableBias = (SqrtTableSize - 1) / max_cell_size2;
		for (int i = 0; i < SqrtTableSize; i++)
		{
			this->SqrtTable[i] = (float)sqrt(i / this->SqrtTableBias);
		}

		this->InputAnalyzedTime.Modified();
	}

	if (renderer->GetRenderWindow()->CheckAbortStatus() || this->GaveError)
	{
		vtkOpenGLCheckErrorMacro("failed during Render");
		return;
	}

	if (renderer->GetRenderWindow()->CheckAbortStatus())
	{
		vtkOpenGLCheckErrorMacro("failed during Render");
		return;
	}

	// Check to see if we need to remap colors.
	if ((this->ColorsMappedTime < this->MTime)
		|| (this->ColorsMappedTime < input->GetMTime())
		|| (this->LastProperty != property)
		|| (this->ColorsMappedTime < property->GetMTime()))
	{
		vtkDataArray *scalars = this->GetScalars(input, this->ScalarMode,
			this->ArrayAccessMode,
			this->ArrayId, this->ArrayName,
			this->UsingCellColors);
		if (!scalars)
		{
			vtkErrorMacro(<< "Can't use projected tetrahedra without scalars!");
			vtkOpenGLCheckErrorMacro("failed during Render");
			return;
		}

		vtkPProjectedTetrahedraMapper::MapScalarsToColors(this->Colors, property,
			scalars);

		this->ColorsMappedTime.Modified();
		this->LastProperty = property;
	}
	if (renderer->GetRenderWindow()->CheckAbortStatus())
	{
		vtkOpenGLCheckErrorMacro("failed during Render");
		return;
	}

	//this->Timer->StartTimer();

	this->ProjectTetrahedra(renderer, volume, renWin);

	//this->Timer->StopTimer();
	//this->TimeToDraw = this->Timer->GetElapsedTime();
	//cout << this->TimeToDraw << endl;
	vtkOpenGLCheckErrorMacro("failed after Render");
}

//-----------------------------------------------------------------------------

inline float vtkOpenGLPProjectedTetrahedraMapper::GetCorrectedDepth(
	float x, float y, float z1, float z2,
	const float inverse_projection_mat[16],
	int use_linear_depth_correction,
	float linear_depth_correction)
{
	if (use_linear_depth_correction)
	{
		float depth = linear_depth_correction*(z1 - z2);
		if (depth < 0) depth = -depth;
		return depth;
	}
	else
	{
		float eye1[3], eye2[3], invw;

		// This code does the same as the commented code above, but also collects
		// common arithmetic between the two matrix x vector operations.  An
		// optimizing compiler may or may not pick up on that.
		float common[4];

		common[0] = (inverse_projection_mat[0] * x
			+ inverse_projection_mat[4] * y
			+ inverse_projection_mat[12]);
		common[1] = (inverse_projection_mat[1] * x
			+ inverse_projection_mat[5] * y
			+ inverse_projection_mat[13]);
		common[2] = (inverse_projection_mat[2] * x
			+ inverse_projection_mat[6] * y
			+ inverse_projection_mat[10] * z1
			+ inverse_projection_mat[14]);
		common[3] = (inverse_projection_mat[3] * x
			+ inverse_projection_mat[7] * y
			+ inverse_projection_mat[15]);

		invw = 1 / (common[3] + inverse_projection_mat[11] * z1);
		eye1[0] = invw*(common[0] + inverse_projection_mat[8] * z1);
		eye1[1] = invw*(common[1] + inverse_projection_mat[9] * z1);
		eye1[2] = invw*(common[2] + inverse_projection_mat[10] * z1);

		invw = 1 / (common[3] + inverse_projection_mat[11] * z2);
		eye2[0] = invw*(common[0] + inverse_projection_mat[8] * z2);
		eye2[1] = invw*(common[1] + inverse_projection_mat[9] * z2);
		eye2[2] = invw*(common[2] + inverse_projection_mat[10] * z2);

		float dist2 = vtkMath::Distance2BetweenPoints(eye1, eye2);
		return this->SqrtTable[(int)(dist2*this->SqrtTableBias)];
	}
}

//-----------------------------------------------------------------------------
void vtkOpenGLPProjectedTetrahedraMapper::ProjectTetrahedra(vtkRenderer *renderer,
	vtkVolume *volume, vtkOpenGLRenderWindow*)
{
	vtkOpenGLClearErrorMacro();

	// after mucking about with FBO bindings be sure
	// we're saving the default fbo attributes/blend function
	this->AllocateFOResources(renderer);

	vtkOpenGLFramebufferObject *fo = NULL;

	// Copy existing Depth/Color  buffers to FO
	if (this->UseFloatingPointFrameBuffer
		&& this->CanDoFloatingPointFrameBuffer)
	{
		fo = this->Framebuffer;

		// bind draw+read to set it up
		fo->SaveCurrentBindingsAndBuffers();
		fo->Bind(fo->GetDrawMode());
		fo->ActivateDrawBuffer(0);

		if (!fo->CheckFrameBufferStatus(fo->GetDrawMode()))
		{
			vtkErrorMacro("FO is incomplete ");
		}

		glBlitFramebuffer(0, 0,
			this->CurrentFBOWidth, this->CurrentFBOHeight,
			0, 0,
			this->CurrentFBOWidth, this->CurrentFBOHeight,
			GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT,
			GL_NEAREST);

		vtkOpenGLCheckErrorMacro("failed at glBlitFramebuffer");
	}

	// TODO:
	// There are some caching optimizations that could be used
	// here to skip various expensive operations (eg sorting
	// cells could be skipped if input data and MVP matrices
	// haven't changed).

	vtkUnstructuredGridBase *input = this->GetInput();

	this->VisibilitySort->SetInput(input);
	this->VisibilitySort->SetDirectionToBackToFront();
	this->VisibilitySort->SetModelTransform(volume->GetMatrix());
	this->VisibilitySort->SetCamera(renderer->GetActiveCamera());
	this->VisibilitySort->SetMaxCellsReturned(1000);
	this->VisibilitySort->InitTraversal();

	this->VisibilitySort->SetNumberOfThread(this->VST);

	this->VisibilitySort->GetCells();
	std::vector<vtkIdTypeArray*> sorted = this->VisibilitySort->GetIDVector();
	std::vector<int> positions = this->VisibilitySort->GetPOSVector();
	if (renderer->GetRenderWindow()->CheckAbortStatus())
	{
		if (fo)
		{
			fo->RestorePreviousBindingsAndBuffers();
		}
		return;
	}
	vtkMatrix4x4 *wcdc;
	vtkMatrix4x4 *wcvc;
	vtkMatrix3x3 *norms;
	vtkMatrix4x4 *vcdc;
	vtkOpenGLCamera *cam = (vtkOpenGLCamera *)(renderer->GetActiveCamera());
	cam->GetKeyMatrices(renderer, wcvc, norms, vcdc, wcdc);
	float projection_mat[16];
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			projection_mat[i * 4 + j] = vcdc->GetElement(i, j);
		}
	}

	float modelview_mat[16];
	if (!volume->GetIsIdentity())
	{
		vtkMatrix4x4 *tmpMat = vtkMatrix4x4::New();
		vtkMatrix4x4 *tmpMat2 = vtkMatrix4x4::New();
		vtkMatrix4x4 *mcwc = volume->GetMatrix();
		tmpMat2->DeepCopy(wcvc);
		tmpMat2->Transpose();
		vtkMatrix4x4::Multiply4x4(tmpMat2, mcwc, tmpMat);
		tmpMat->Transpose();
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				modelview_mat[i * 4 + j] = tmpMat->GetElement(i, j);
			}
		}
		tmpMat->Delete();
		tmpMat2->Delete();
	}
	else
	{
		for (int i = 0; i < 4; ++i)
		{
			for (int j = 0; j < 4; ++j)
			{
				modelview_mat[i * 4 + j] = wcvc->GetElement(i, j);
			}
		}
	}

	// Get the inverse projection matrix so that we can convert distances in
	// clipping space to distances in world or eye space.
	float inverse_projection_mat[16];
	float linear_depth_correction = 1;
	int use_linear_depth_correction;
	double tmp_mat[16];

	// VTK's matrix functions use doubles.
	std::copy(projection_mat, projection_mat + 16, tmp_mat);
	// VTK and OpenGL store their matrices differently.  Correct.
	vtkMatrix4x4::Transpose(tmp_mat, tmp_mat);
	// Take the inverse.
	vtkMatrix4x4::Invert(tmp_mat, tmp_mat);
	// Restore back to OpenGL form.
	vtkMatrix4x4::Transpose(tmp_mat, tmp_mat);
	// Copy back to float for faster computation.
	std::copy(tmp_mat, tmp_mat + 16, inverse_projection_mat);

	// Check to see if we can just do a linear depth correction from clipping
	// space to eye space.
	use_linear_depth_correction = ((projection_mat[3] == 0.0)
		&& (projection_mat[7] == 0.0)
		&& (projection_mat[11] == 0.0)
		&& (projection_mat[15] == 1.0));
	if (use_linear_depth_correction)
	{
		float pos1[3], *pos2;

		pos1[0] = inverse_projection_mat[8] + inverse_projection_mat[12];
		pos1[1] = inverse_projection_mat[9] + inverse_projection_mat[13];
		pos1[2] = inverse_projection_mat[10] + inverse_projection_mat[14];

		pos2 = inverse_projection_mat + 12;

		linear_depth_correction = sqrt(vtkMath::Distance2BetweenPoints(pos1, pos2));
	}
	// Transform all the points.
	vtkPProjectedTetrahedraMapper::TransformPoints(input->GetPoints(),
		projection_mat, modelview_mat,
		this->TransformedPoints);
	float *points = this->TransformedPoints->GetPointer(0);

	if (renderer->GetRenderWindow()->CheckAbortStatus())
	{
		if (fo)
		{
			fo->RestorePreviousBindingsAndBuffers();
		}
		return;
	}

	glDepthMask(GL_FALSE);

	glDisable(GL_CULL_FACE);

	GLint blendSrcA = GL_ONE;
	GLint blendDstA = GL_ONE_MINUS_SRC_ALPHA;
	GLint blendSrcC = GL_SRC_ALPHA;
	GLint blendDstC = GL_ONE_MINUS_SRC_ALPHA;
	glGetIntegerv(GL_BLEND_SRC_ALPHA, &blendSrcA);
	glGetIntegerv(GL_BLEND_DST_ALPHA, &blendDstA);
	glGetIntegerv(GL_BLEND_SRC_RGB, &blendSrcC);
	glGetIntegerv(GL_BLEND_DST_RGB, &blendDstC);
	glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA,
		GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

	float unit_distance = volume->GetProperty()->GetScalarOpacityUnitDistance();

	// build the VBO and IBOs,  we so these in chuncks as based on
	// the settings of the VisibilitySort tclass
	this->VBO->SetStride(6 * sizeof(float));


	unsigned char *colors = this->Colors->GetPointer(0);
	vtkIdType totalnumcells = input->GetNumberOfCells();
	vtkIdType numcellsrendered = 0;

	//std::vector<float> packedVBO;
	//packedVBO.reserve(6 * 5 * this->VisibilitySort->GetMaxCellsReturned());

	//std::vector<unsigned int> indexArray;
	//indexArray.reserve(3 * 4 * this->VisibilitySort->GetMaxCellsReturned());
	// Let's do it!
#pragma region kernel loop
	//for (vtkIdTypeArray *sorted_cell_ids = this->VisibilitySort->GetNextCells();
	//	sorted_cell_ids != NULL;
	//	sorted_cell_ids = this->VisibilitySort->GetNextCells())

	this->Timer->StartTimer();
#pragma omp parallel for num_threads(this->CDT)//设置线程数
	for (int s = 0; s < sorted.size(); ++s)
	{
		/*this->UpdateProgress((double)numcellsrendered / totalnumcells);
		if (renderer->GetRenderWindow()->CheckAbortStatus())
		{
		break;
		}*/

		//vtkIdType *cell_ids = sorted_cell_ids->GetPointer(0);
		//vtkIdType num_cell_ids = sorted_cell_ids->GetNumberOfTuples();

		vtkIdType *cell_ids = sorted[s]->GetPointer(0);
		vtkIdType num_cell_ids = sorted[s]->GetNumberOfTuples();

		unsigned int *currentII = this->indexArray + 12 * positions[s];
		float *currentPI = this->packedVBO + 30 * positions[s];
		/*	packedVBO.resize(6 * 5 * num_cell_ids);
		std::vector<float>::iterator it = packedVBO.begin();*/

		int numPts = 0;
		/*indexArray.resize(0);*/
		/*	vtkNew<vtkIdList> cellPointIds;*/
#pragma region single order visibility cells
		//#pragma omp parallel for private(cellPointIds)  num_threads(2)//设置线程数
		//#pragma omp parallel for num_threads(20)//设置线程数
		for (vtkIdType i = 0; i < num_cell_ids; i++)
		{
			float tet_points[5 * 3] = { 0.0f };
			unsigned char tet_colors[5 * 3] = { 0 };
			float tet_texcoords[5 * 2] = { 0.0f };

			vtkIdType cell = cell_ids[i];
			//input->GetCellPoints(cell, cellPointIds.GetPointer());
			unsigned int * cellPoints = this->cellArray + 4 * cell;

			int j;

			// Get the data for the tetrahedra.
			for (j = 0; j < 4; j++)
			{
				// Assuming we only have tetrahedra, each entry in cells has 5
				// components.
				//const float *p = points + 3 * cellPointIds->GetId(j);
				const float *p = points + 3 * cellPoints[j];
				tet_points[j * 3 + 0] = p[0];
				tet_points[j * 3 + 1] = p[1];
				tet_points[j * 3 + 2] = p[2];

				const unsigned char *c;
				if (this->UsingCellColors)
				{
					c = colors + 4 * cell;
				}
				else
				{
					//c = colors + 4 * cellPointIds->GetId(j);
					c = colors + 4 * cellPoints[j];
				}

				tet_colors[j * 3 + 0] = c[0];
				tet_colors[j * 3 + 1] = c[1];
				tet_colors[j * 3 + 2] = c[2];

				tet_texcoords[j * 2 + 0] = static_cast<float>(c[3]) / 255.0f;
				tet_texcoords[j * 2 + 1] = 0;
			}

			unsigned int *tmpII = currentII + 12 * i;
			float *tmpPI = currentPI + 30 * i;

			if (((tet_points[0 * 3 + 0] > 1.0f) && (tet_points[1 * 3 + 0] > 1.0f)
				&& (tet_points[2 * 3 + 0] > 1.0f) && (tet_points[3 * 3 + 0] > 1.0f))
				|| ((tet_points[0 * 3 + 0] < -1.0f) && (tet_points[1 * 3 + 0] < -1.0f)
					&& (tet_points[2 * 3 + 0] < -1.0f) && (tet_points[3 * 3 + 0] < -1.0f))
				|| ((tet_points[0 * 3 + 1] > 1.0f) && (tet_points[1 * 3 + 1] > 1.0f)
					&& (tet_points[2 * 3 + 1] > 1.0f) && (tet_points[3 * 3 + 1] > 1.0f))
				|| ((tet_points[0 * 3 + 1] < -1.0f) && (tet_points[1 * 3 + 1] < -1.0f)
					&& (tet_points[2 * 3 + 1] < -1.0f) && (tet_points[3 * 3 + 1] < -1.0f))
				|| ((tet_points[0 * 3 + 2] > 1.0f) && (tet_points[1 * 3 + 2] > 1.0f)
					&& (tet_points[2 * 3 + 2] > 1.0f) && (tet_points[3 * 3 + 2] > 1.0f))
				|| ((tet_points[0 * 3 + 2] < -1.0f) || (tet_points[1 * 3 + 2] < -1.0f)
					|| (tet_points[2 * 3 + 2] < -1.0f) || (tet_points[3 * 3 + 2] < -1.0f)))
			{
				for (int cellIdx = 0; cellIdx < 4; cellIdx++)
				{
					*(tmpII++) = *(tmpII++) = *(tmpII++) = 5 * i + 5 * positions[s];
				}
				for (int ptIdx = 0; ptIdx < 5; ptIdx++)
				{
					*(tmpPI++) = tet_points[ptIdx * 3];
					*(tmpPI++) = tet_points[ptIdx * 3 + 1];
					*(tmpPI++) = tet_points[ptIdx * 3 + 2];
					*(tmpPI++) = 0;
					*(tmpPI++) = 0; // attenuation
					*(tmpPI++) = 0; // depth
				}
			}
			else
			{
				vtkIdType segment1[2];
				vtkIdType segment2[2];
				float v1[2], v2[2], v3[3];
				v1[0] = tet_points[1 * 3 + 0] - tet_points[0 * 3 + 0];
				v1[1] = tet_points[1 * 3 + 1] - tet_points[0 * 3 + 1];
				v2[0] = tet_points[2 * 3 + 0] - tet_points[0 * 3 + 0];
				v2[1] = tet_points[2 * 3 + 1] - tet_points[0 * 3 + 1];
				v3[0] = tet_points[3 * 3 + 0] - tet_points[0 * 3 + 0];
				v3[1] = tet_points[3 * 3 + 1] - tet_points[0 * 3 + 1];

				float face_dir1 = v3[0] * v2[1] - v3[1] * v2[0];
				float face_dir2 = v1[0] * v3[1] - v1[1] * v3[0];
				float face_dir3 = v2[0] * v1[1] - v2[1] * v1[0];

				if ((face_dir1 * face_dir2 >= 0)
					&& ((face_dir1 != 0)       // Handle a special case where 2 faces
						|| (face_dir2 != 0)))   // are perpendicular to the view plane.
				{
					segment1[0] = 0;  segment1[1] = 3;
					segment2[0] = 1;  segment2[1] = 2;
				}
				else if (face_dir1 * face_dir3 >= 0)
				{
					segment1[0] = 0;  segment1[1] = 2;
					segment2[0] = 1;  segment2[1] = 3;
				}
				else      // Unless the tet is degenerate, face_dir2*face_dir3 >= 0
				{
					segment1[0] = 0;  segment1[1] = 1;
					segment2[0] = 2;  segment2[1] = 3;
				}

#define VEC3SUB(Z,X,Y)          \
  (Z)[0] = (X)[0] - (Y)[0];     \
  (Z)[1] = (X)[1] - (Y)[1];     \
  (Z)[2] = (X)[2] - (Y)[2];
#define P1 (tet_points + 3*segment1[0])
#define P2 (tet_points + 3*segment1[1])
#define P3 (tet_points + 3*segment2[0])
#define P4 (tet_points + 3*segment2[1])
#define C1 (tet_colors + 3*segment1[0])
#define C2 (tet_colors + 3*segment1[1])
#define C3 (tet_colors + 3*segment2[0])
#define C4 (tet_colors + 3*segment2[1])
#define T1 (tet_texcoords + 2*segment1[0])
#define T2 (tet_texcoords + 2*segment1[1])
#define T3 (tet_texcoords + 2*segment2[0])
#define T4 (tet_texcoords + 2*segment2[1])
				float A[3], B[3], C[3];
				VEC3SUB(A, P2, P1);
				VEC3SUB(B, P4, P3);
				VEC3SUB(C, P3, P1);
				float denominator = (A[0] * B[1] - A[1] * B[0]);
				if (denominator == 0) continue;   // Must be degenerated tetrahedra.
				float alpha = (B[1] * C[0] - B[0] * C[1]) / denominator;
				float beta = (A[1] * C[0] - A[0] * C[1]) / denominator;

				if ((alpha >= 0) && (alpha <= 1))
				{
					tet_points[3 * 4 + 0] = P1[0] + alpha*A[0];
					tet_points[3 * 4 + 1] = P1[1] + alpha*A[1];
					tet_points[3 * 4 + 2] = P1[2] + alpha*A[2];

					float depth = this->GetCorrectedDepth(
						tet_points[3 * 4 + 0],
						tet_points[3 * 4 + 1],
						tet_points[3 * 4 + 2],
						P3[2] + beta*B[2],
						inverse_projection_mat,
						use_linear_depth_correction,
						linear_depth_correction);

					tet_colors[3 * 4 + 0] = static_cast<unsigned char>
						(0.5f*(C1[0] + alpha*(C2[0] - C1[0])
							+ C3[0] + beta*(C4[0] - C3[0])));

					tet_colors[3 * 4 + 1] = static_cast<unsigned char>
						(0.5f*(C1[1] + alpha*(C2[1] - C1[1])
							+ C3[1] + beta*(C4[1] - C3[1])));

					tet_colors[3 * 4 + 2] = static_cast<unsigned char>
						(0.5f*(C1[2] + alpha*(C2[2] - C1[2])
							+ C3[2] + beta*(C4[2] - C3[2])));

					tet_texcoords[2 * 4 + 0] = 0.5f*(T1[0] + alpha*(T2[0] - T1[0])
						+ T3[0] + alpha*(T4[0] - T3[0]));

					tet_texcoords[2 * 4 + 1] = depth / unit_distance;

					unsigned char indices[6];
					indices[0] = 4;
					indices[1] = segment1[0];
					indices[2] = segment2[0];
					indices[3] = segment1[1];
					indices[4] = segment2[1];
					indices[5] = segment1[0];
					for (int cellIdx = 0; cellIdx < 4; cellIdx++)
					{
						*(tmpII++) = indices[0] + 5 * i + 5 * positions[s];
						*(tmpII++) = indices[cellIdx + 1] + 5 * i + 5 * positions[s];
						*(tmpII++) = indices[cellIdx + 2] + 5 * i + 5 * positions[s];


						/*indexArray.push_back(indices[0] + numPts);
						indexArray.push_back(indices[cellIdx + 1] + numPts);
						indexArray.push_back(indices[cellIdx + 2] + numPts);*/
					}
				}
				else
				{
					if (alpha <= 0)
					{
						std::swap(segment1[0], segment1[1]);
						alpha = 1 - alpha;
					}
					float edgez = P3[2] + beta*B[2];
					float pointz = P1[2];
					float facez = (edgez + (alpha - 1)*pointz) / alpha;
					float depth = GetCorrectedDepth(P2[0], P2[1], P2[2], facez,
						inverse_projection_mat,
						use_linear_depth_correction,
						linear_depth_correction);
					for (j = 0; j < 3; j++)
					{
						float edgec = C3[j] + beta*(C4[j] - C3[j]);
						float pointc = C1[j];
						float facec = (edgec + (alpha - 1)*pointc) / alpha;
						C2[j] = (unsigned char)(0.5f*(facec + C2[j]));
					}
					float edgea = T3[0] + beta*(T4[0] - T3[0]);
					float pointa = T1[0];
					float facea = (edgea + (alpha - 1)*pointa) / alpha;
					T2[0] = 0.5f*(facea + T2[0]);
					T2[1] = depth / unit_distance;

					unsigned char indices[5];
					indices[0] = segment1[1];
					indices[1] = segment1[0];
					indices[2] = segment2[0];
					indices[3] = segment2[1];
					indices[4] = segment1[0];

					*(tmpII++) = indices[0] + 5 * i + 5 * positions[s];
					*(tmpII++) = indices[1] + 5 * i + 5 * positions[s];
					*(tmpII++) = indices[2] + 5 * i + 5 * positions[s];

					*(tmpII++) = indices[0] + 5 * i + 5 * positions[s];
					*(tmpII++) = indices[2] + 5 * i + 5 * positions[s];
					*(tmpII++) = indices[3] + 5 * i + 5 * positions[s];

					*(tmpII++) = indices[0] + 5 * i + 5 * positions[s];
					*(tmpII++) = indices[3] + 5 * i + 5 * positions[s];
					*(tmpII++) = indices[4] + 5 * i + 5 * positions[s];

					//for (int cellIdx = 0; cellIdx < 3; cellIdx++)
					//{
					//	*(tmpII++) = indices[0] + 5 * i + 5 * positions[s];
					//	*(tmpII++) = indices[cellIdx + 1] + 5 * i + 5 * positions[s];
					//	*(tmpII++) = indices[cellIdx + 2] + 5 * i + 5 * positions[s];

					//	/*	indexArray.push_back(indices[0] + numPts);
					//	indexArray.push_back(indices[cellIdx + 1] + numPts);
					//	indexArray.push_back(indices[cellIdx + 2] + numPts);*/
					//}
					*(tmpII++) = indices[0] + 5 * i + 5 * positions[s];
					*(tmpII++) = indices[0] + 5 * i + 5 * positions[s];
					*(tmpII++) = indices[0] + 5 * i + 5 * positions[s];

				}

				union { unsigned char c[4]; float f; } v = { { 0, 0, 0, 255 } };
				for (int ptIdx = 0; ptIdx < 5; ptIdx++)
				{
					*(tmpPI++) = tet_points[ptIdx * 3];
					*(tmpPI++) = tet_points[ptIdx * 3 + 1];
					*(tmpPI++) = tet_points[ptIdx * 3 + 2];
					v.c[0] = tet_colors[ptIdx * 3];
					v.c[1] = tet_colors[ptIdx * 3 + 1];
					v.c[2] = tet_colors[ptIdx * 3 + 2];
					*(tmpPI++) = v.f;
					*(tmpPI++) = tet_texcoords[ptIdx * 2]; // attenuation
					*(tmpPI++) = tet_texcoords[ptIdx * 2 + 1]; // depth
				}
			}

			//numPts = 5 * (i + 1);
		}
#pragma endregion
	}
#pragma endregion

	this->Timer->StopTimer();
	cout << "classification:\t" << this->Timer->GetElapsedTime() << endl;

	this->Timer->StartTimer();
	this->VBO->Upload(packedVBO, 30 * input->GetNumberOfCells(), vtkOpenGLBufferObject::ArrayBuffer);
	this->VBO->Bind();

	this->Tris.VAO->Bind();
	if (this->Tris.IBO->IndexCount && (
		this->Tris.ShaderSourceTime > this->Tris.AttributeUpdateTime))
	{
		if (!this->Tris.VAO->AddAttributeArray(this->Tris.Program, this->VBO,
			"vertexDC", 0,
			this->VBO->GetStride(), VTK_FLOAT, 3, false))
		{
			vtkErrorMacro(<< "Error setting 'vertexDC' in shader VAO.");
		}
		if (!this->Tris.VAO->AddAttributeArray(this->Tris.Program, this->VBO,
			"scalarColor", 3 * sizeof(float),
			this->VBO->GetStride(), VTK_UNSIGNED_CHAR,
			3, true))
		{
			vtkErrorMacro(<< "Error setting 'scalarColor' in shader VAO.");
		}
		if (!this->Tris.VAO->AddAttributeArray(this->Tris.Program, this->VBO,
			"attenuationArray", 4 * sizeof(float),
			this->VBO->GetStride(), VTK_FLOAT,
			1, false))
		{
			vtkErrorMacro(<< "Error setting attenuation in shader VAO.");
		}
		if (!this->Tris.VAO->AddAttributeArray(this->Tris.Program, this->VBO,
			"depthArray", 5 * sizeof(float),
			this->VBO->GetStride(), VTK_FLOAT,
			1, false))
		{
			vtkErrorMacro(<< "Error setting depth in shader VAO.");
		}
		this->Tris.AttributeUpdateTime.Modified();
	}


	//this->Tris.IBO->Upload(indexArray, vtkOpenGLBufferObject::ElementArrayBuffer);
	//this->Tris.IBO->IndexCount = indexArray.size();

	this->Tris.IBO->Upload(indexArray, 12 * input->GetNumberOfCells(), vtkOpenGLBufferObject::ElementArrayBuffer);
	this->Tris.IBO->IndexCount = 12 * input->GetNumberOfCells();

	this->Tris.IBO->Bind();
	glDrawRangeElements(GL_TRIANGLES, 0,
		static_cast<GLuint>(12 * input->GetNumberOfCells() - 1),
		static_cast<GLsizei>(this->Tris.IBO->IndexCount),
		GL_UNSIGNED_INT,
		reinterpret_cast<const GLvoid *>(NULL));
	this->Tris.IBO->Release();
	this->Tris.VAO->Release();
	this->VBO->Release();
	//numcellsrendered += num_cell_ids;

	if (fo)
	{
		// copy from our fbo to the default one
		fo->Bind(fo->GetReadMode());

		// draw to default fbo
		fo->RestorePreviousBindingsAndBuffers(fo->GetDrawMode());

		// Depth buffer has not changed so only copy color
		glBlitFramebuffer(0, 0, this->CurrentFBOWidth, this->CurrentFBOHeight,
			0, 0, this->CurrentFBOWidth, this->CurrentFBOHeight,
			GL_COLOR_BUFFER_BIT, GL_NEAREST);

		vtkOpenGLCheckErrorMacro("failed at glBlitFramebuffer");

		// restore default fbo for both read+draw
		fo->RestorePreviousBindingsAndBuffers(fo->GetReadMode());
	}

	// Restore the blend function.
	vtkOpenGLCheckErrorMacro("failed at glPopAttrib");

	glDepthMask(GL_TRUE);
	glBlendFuncSeparate(blendSrcC, blendDstC, blendSrcA, blendDstA);

	vtkOpenGLCheckErrorMacro("failed after ProjectTetrahedra");
	this->UpdateProgress(1.0);

	this->Timer->StopTimer();
	cout << "raster:\t" << this->Timer->GetElapsedTime() << endl;
}