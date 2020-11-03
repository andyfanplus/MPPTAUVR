#ifndef vtkPProjectedTetrahedraMapper_h
#define vtkPProjectedTetrahedraMapper_h

#include "vtkRenderingVolumeModule.h" // For export macro
#include "vtkUnstructuredGridVolumeMapper.h"

class vtkFloatArray;
class vtkPoints;
class vtkUnsignedCharArray;
class vtkPVisibilitySort;//
class vtkVolumeProperty;
class vtkRenderWindow;

class  vtkPProjectedTetrahedraMapper : public vtkUnstructuredGridVolumeMapper
{
public:
  virtual	void SetNumberOfVST(int)=0;
  virtual	void SetNumberOfCDT(int) = 0;

  vtkTypeMacro(vtkPProjectedTetrahedraMapper,
                       vtkUnstructuredGridVolumeMapper);
  static vtkPProjectedTetrahedraMapper *New();
  void PrintSelf(ostream &os, vtkIndent indent) VTK_OVERRIDE;

  virtual void SetVisibilitySort(vtkPVisibilitySort *sort);
  vtkGetObjectMacro(VisibilitySort, vtkPVisibilitySort);

  static void MapScalarsToColors(vtkDataArray *colors,
                                 vtkVolumeProperty *property,
                                 vtkDataArray *scalars);
  static void TransformPoints(vtkPoints *inPoints,
                              const float projection_mat[16],
                              const float modelview_mat[16],
                              vtkFloatArray *outPoints);

  /**
   * Return true if the rendering context provides
   * the nececessary functionality to use this class.
   */
  virtual bool IsSupported(vtkRenderWindow *)
    { return false; }


protected:
  vtkPProjectedTetrahedraMapper();
  ~vtkPProjectedTetrahedraMapper() VTK_OVERRIDE;

  vtkPVisibilitySort *VisibilitySort;

  /**
   * The visibility sort will probably make a reference loop by holding a
   * reference to the input.
   */
  void ReportReferences(vtkGarbageCollector *collector) VTK_OVERRIDE;

  int VST = 1;
  int CDT = 1;

private:
  vtkPProjectedTetrahedraMapper(const vtkPProjectedTetrahedraMapper &) VTK_DELETE_FUNCTION;
  void operator=(const vtkPProjectedTetrahedraMapper &) VTK_DELETE_FUNCTION;
};

#endif
