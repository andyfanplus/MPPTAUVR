
#ifndef vtkOpenGLPProjectedTetrahedraMapper_h
#define vtkOpenGLPProjectedTetrahedraMapper_h

#include "vtkRenderingVolumeOpenGL2Module.h" // For export macro
#include "vtkPProjectedTetrahedraMapper.h"

#include "vtkOpenGLHelper.h" // used for ivars

class vtkPVisibilitySort;
class vtkUnsignedCharArray;
class vtkFloatArray;
class vtkRenderWindow;
class vtkOpenGLFramebufferObject;
class vtkOpenGLRenderWindow;
class vtkOpenGLVertexBufferObject;

class  vtkOpenGLPProjectedTetrahedraMapper
  : public vtkPProjectedTetrahedraMapper
{
public:
  vtkTypeMacro(vtkOpenGLPProjectedTetrahedraMapper,
                       vtkPProjectedTetrahedraMapper);
  static vtkOpenGLPProjectedTetrahedraMapper *New();
  void PrintSelf(ostream &os, vtkIndent indent) VTK_OVERRIDE;

  void ReleaseGraphicsResources(vtkWindow *window) VTK_OVERRIDE;

  void Render(vtkRenderer *renderer, vtkVolume *volume) VTK_OVERRIDE;

  //@{
  /**
   * Set/get whether to use floating-point rendering buffers rather
   * than the default.
   */
  vtkSetMacro(UseFloatingPointFrameBuffer,bool);
  vtkGetMacro(UseFloatingPointFrameBuffer,bool);
  vtkBooleanMacro(UseFloatingPointFrameBuffer,bool);
  //@}

  /**
   * Return true if the rendering context provides
   * the nececessary functionality to use this class.
   */
  bool IsSupported(vtkRenderWindow *context) VTK_OVERRIDE;

  virtual void  SetNumberOfVST(int vst);
  virtual void  SetNumberOfCDT(int vst);
protected:
  vtkOpenGLPProjectedTetrahedraMapper();
  ~vtkOpenGLPProjectedTetrahedraMapper() VTK_OVERRIDE;

  void Initialize(vtkRenderer *ren);
  bool Initialized;
  int  CurrentFBOWidth, CurrentFBOHeight;
  bool AllocateFOResources(vtkRenderer *ren);
  bool CanDoFloatingPointFrameBuffer;
  bool FloatingPointFrameBufferResourcesAllocated;
  bool UseFloatingPointFrameBuffer;
  bool HasHardwareSupport;

  vtkUnsignedCharArray *Colors;
  int UsingCellColors;

  vtkFloatArray *TransformedPoints;

  float MaxCellSize;
  vtkTimeStamp InputAnalyzedTime;
  vtkTimeStamp ColorsMappedTime;

  // The VBO and its layout.
  vtkOpenGLVertexBufferObject *VBO;

  // Structures for the various cell types we render.
  vtkOpenGLHelper Tris;

  int GaveError;

  vtkVolumeProperty *LastProperty;

  vtkOpenGLFramebufferObject *Framebuffer;

  float *SqrtTable;
  float SqrtTableBias;

  virtual void ProjectTetrahedra(vtkRenderer *renderer, vtkVolume *volume,
    vtkOpenGLRenderWindow* renWin);

  float GetCorrectedDepth(float x, float y, float z1, float z2,
                          const float inverse_projection_mat[16],
                          int use_linear_depth_correction,
                          float linear_depth_correction);

private:
  vtkOpenGLPProjectedTetrahedraMapper(const vtkOpenGLPProjectedTetrahedraMapper &) VTK_DELETE_FUNCTION;
  void operator=(const vtkOpenGLPProjectedTetrahedraMapper &) VTK_DELETE_FUNCTION;

  class vtkInternals;
  vtkInternals *Internals;

  float *packedVBO;
  unsigned int *indexArray;
  unsigned int *cellArray;

};

#endif
