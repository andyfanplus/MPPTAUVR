#ifndef vtkPVisibilitySort_h
#define vtkPVisibilitySort_h

#include "vtkRenderingCoreModule.h" // For export macro
#include "vtkObject.h"
#include "vector"

class vtkIdTypeArray;
class vtkDataSet;
class vtkMatrix4x4;
class vtkCamera;
class vector;

using namespace std;
class  vtkPVisibilitySort : public vtkObject
{
public:
  vtkTypeMacro(vtkPVisibilitySort, vtkObject);
  void PrintSelf(ostream &os, vtkIndent indent) VTK_OVERRIDE;

  //@{
  /**
   * To facilitate incremental sorting algorithms, the cells are retrieved
   * in an iteration process.  That is, call InitTraversal to start the
   * iteration and call GetNextCells to get the cell IDs in order.
   * However, for efficiencies sake, GetNextCells returns an ordered list
   * of several id's in once call (but not necessarily all).  GetNextCells
   * will return NULL once the entire sorted list is output.  The
   * vtkIdTypeArray returned from GetNextCells is a cached array, so do not
   * delete it.  At the same note, do not expect the array to be valid
   * after subsequent calls to GetNextCells.
   */
  virtual void InitTraversal() = 0;
  virtual vtkIdTypeArray *GetNextCells() = 0;
  //@}

  //@{
  /**
   * Set/Get the maximum number of cells that GetNextCells will return
   * in one invocation.
   */
  vtkSetClampMacro(MaxCellsReturned, int, 1, VTK_INT_MAX);
  vtkGetMacro(MaxCellsReturned, int);
  //@}

  //@{
  /**
   * Set/Get the matrix that transforms from object space to world space.
   * Generally, you get this matrix from a call to GetMatrix of a vtkProp3D
   * (vtkActor).
   */
  virtual void SetModelTransform(vtkMatrix4x4 *mat);
  vtkGetObjectMacro(ModelTransform, vtkMatrix4x4);
  //@}

  vtkGetObjectMacro(InverseModelTransform, vtkMatrix4x4);

  //@{
  /**
   * Set/Get the camera that specifies the viewing parameters.
   */
  virtual void SetCamera(vtkCamera *camera);
  vtkGetObjectMacro(Camera, vtkCamera);
  //@}

  //@{
  /**
   * Set/Get the data set containing the cells to sort.
   */
  virtual void SetInput(vtkDataSet *data);
  vtkGetObjectMacro(Input, vtkDataSet);
  //@}

  //@{
  /**
   * Set/Get the sorting direction.  Be default, the direction is set
   * to back to front.
   */
  vtkGetMacro(Direction, int);
  vtkSetMacro(Direction, int);
  void SetDirectionToBackToFront() { this->SetDirection(BACK_TO_FRONT); }
  void SetDirectionToFrontToBack() { this->SetDirection(FRONT_TO_BACK); }
  //@}

  enum { BACK_TO_FRONT, FRONT_TO_BACK };

  //@{
  /**
   * Overwritten to enable garbage collection.
   */
  void Register(vtkObjectBase *o) VTK_OVERRIDE;
  void UnRegister(vtkObjectBase *o) VTK_OVERRIDE;
  //@}

  virtual void GetCells() = 0;
  virtual std::vector<vtkIdTypeArray*> GetIDVector() = 0;
  virtual std::vector<int> GetPOSVector() = 0;

  virtual void SetNumberOfThread(int) = 0;
protected:
  vtkPVisibilitySort();
  ~vtkPVisibilitySort() VTK_OVERRIDE;

  vtkTimeStamp LastSortTime;

  vtkMatrix4x4 *ModelTransform;
  vtkMatrix4x4 *InverseModelTransform;
  vtkCamera *Camera;
  vtkDataSet *Input;

  int MaxCellsReturned;

  int Direction;

  void ReportReferences(vtkGarbageCollector *collector) VTK_OVERRIDE;

  int Thread=1;

private:
  vtkPVisibilitySort(const vtkPVisibilitySort &) VTK_DELETE_FUNCTION;
  void operator=(const vtkPVisibilitySort &) VTK_DELETE_FUNCTION;
};

#endif //vtkPVisibilitySort_h

