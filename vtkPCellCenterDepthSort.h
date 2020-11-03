#ifndef vtkPCellCenterDepthSort_h
#define vtkPCellCenterDepthSort_h

#include "vtkRenderingCoreModule.h" // For export macro
#include "vtkPVisibilitySort.h"
class vtkFloatArray;
class vtkPCellCenterDepthSortStack;
class vtkTimerLog;

class  vtkPCellCenterDepthSort : public vtkPVisibilitySort
{
public:
  vtkTypeMacro(vtkPCellCenterDepthSort, vtkPVisibilitySort);
  void PrintSelf(ostream &os, vtkIndent indent) VTK_OVERRIDE;
  static vtkPCellCenterDepthSort *New();

  void InitTraversal() VTK_OVERRIDE;
  vtkIdTypeArray *GetNextCells() VTK_OVERRIDE;
  void GetCells() VTK_OVERRIDE;
  std::vector<vtkIdTypeArray*> GetIDVector() VTK_OVERRIDE;
  std::vector<int> GetPOSVector()VTK_OVERRIDE;

  virtual void SetNumberOfThread(int);
protected:
  vtkPCellCenterDepthSort();
  ~vtkPCellCenterDepthSort() VTK_OVERRIDE;

  vtkIdTypeArray *SortedCells;
  vtkIdTypeArray *SortedCellPartition;

  vtkFloatArray *CellCenters;
  vtkFloatArray *CellDepths;
  vtkFloatArray *CellPartitionDepths;

  vtkTimerLog *Timer;

  virtual float *ComputeProjectionVector();
  virtual void ComputeCellCenters();
  virtual void ComputeDepths();

private:
  vtkPCellCenterDepthSortStack *ToSort;

  vtkPCellCenterDepthSort(const vtkPCellCenterDepthSort &) VTK_DELETE_FUNCTION;
  void operator=(const vtkPCellCenterDepthSort &) VTK_DELETE_FUNCTION;
};

#endif
