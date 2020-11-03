#include "vtkPVisibilitySort.h"

#include "vtkIdList.h"
#include "vtkDataSet.h"
#include "vtkMatrix4x4.h"
#include "vtkCamera.h"
#include "vtkGarbageCollector.h"

//-----------------------------------------------------------------------------

vtkCxxSetObjectMacro(vtkPVisibilitySort, Camera, vtkCamera);
vtkCxxSetObjectMacro(vtkPVisibilitySort, Input, vtkDataSet);

//-----------------------------------------------------------------------------

vtkPVisibilitySort::vtkPVisibilitySort()
{
  this->ModelTransform = vtkMatrix4x4::New();
  this->ModelTransform->Identity();
  this->InverseModelTransform = vtkMatrix4x4::New();
  this->InverseModelTransform->Identity();

  this->Camera = NULL;
  this->Input = NULL;

  this->Direction = vtkPVisibilitySort::BACK_TO_FRONT;

  this->MaxCellsReturned = VTK_INT_MAX;
}

//-----------------------------------------------------------------------------

vtkPVisibilitySort::~vtkPVisibilitySort()
{
  this->ModelTransform->Delete();
  this->InverseModelTransform->Delete();

  this->SetCamera(NULL);
  this->SetInput(NULL);
}

//-----------------------------------------------------------------------------

void vtkPVisibilitySort::Register(vtkObjectBase *o)
{
  this->RegisterInternal(o, 1);
}

void vtkPVisibilitySort::UnRegister(vtkObjectBase *o)
{
  this->UnRegisterInternal(o, 1);
}

void vtkPVisibilitySort::ReportReferences(vtkGarbageCollector *collector)
{
  this->Superclass::ReportReferences(collector);
  vtkGarbageCollectorReport(collector, this->Input, "Input");
}

//-----------------------------------------------------------------------------

void vtkPVisibilitySort::SetModelTransform(vtkMatrix4x4 *mat)
{
  // Less efficient than vtkMatrix4x4::DeepCopy, but only sets Modified if
  // there is a real change.
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      this->ModelTransform->SetElement(i, j, mat->GetElement(i, j));
    }
  }

  if (  this->ModelTransform->GetMTime()
      > this->InverseModelTransform->GetMTime() )
  {
    this->InverseModelTransform->DeepCopy(this->ModelTransform);
    this->InverseModelTransform->Invert();
  }
}

//-----------------------------------------------------------------------------

void vtkPVisibilitySort::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << "Input: (" << this->Input << ")" << endl;
  os << indent << "Direction: ";
  switch (this->Direction)
  {
    case vtkPVisibilitySort::BACK_TO_FRONT:
      os << "back to front" << endl;
      break;
    case vtkPVisibilitySort::FRONT_TO_BACK:
      os << "front to back" << endl;
      break;
    default:
      os << "unknown" << endl;
      break;
  }

  os << indent << "MaxCellsReturned: " << this->MaxCellsReturned << endl;

  os << indent << "ModelTransform:" << endl;
  this->ModelTransform->PrintSelf(os, indent.GetNextIndent());
  os << indent << "InverseModelTransform:" << endl;
  this->InverseModelTransform->PrintSelf(os, indent.GetNextIndent());

  os << indent << "Camera: (" << this->Camera << ")" << endl;
}
