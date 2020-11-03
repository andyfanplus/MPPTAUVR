#include "vtkPCellCenterDepthSort.h"

#include "vtkObjectFactory.h"
#include "vtkIdTypeArray.h"
#include "vtkDataSet.h"
#include "vtkCamera.h"
#include "vtkMatrix4x4.h"
#include "vtkFloatArray.h"
#include "vtkCell.h"
#include "vtkMath.h"
#include "vtkSortDataArray.h"

#include"vtkTimerLog.h"

#include <stack>
#include <utility>
#include <algorithm>
#include<omp.h>
//-----------------------------------------------------------------------------

typedef std::pair<vtkIdType, vtkIdType> vtkIdPair;

class vtkPCellCenterDepthSortStack
{
public:
  std::stack<vtkIdPair> Stack;
  std::vector<vtkIdPair> Vector;
  std::vector<vtkIdTypeArray*> IDVector;
  std::vector<int> POSVector;
};

//-----------------------------------------------------------------------------

vtkStandardNewMacro(vtkPCellCenterDepthSort);

vtkPCellCenterDepthSort::vtkPCellCenterDepthSort()
{
  this->SortedCells = vtkIdTypeArray::New();
  this->SortedCells->SetNumberOfComponents(1);
  this->SortedCellPartition = vtkIdTypeArray::New();
  this->SortedCells->SetNumberOfComponents(1);

  this->CellCenters = vtkFloatArray::New();
  this->CellCenters->SetNumberOfComponents(3);
  this->CellDepths = vtkFloatArray::New();
  this->CellDepths->SetNumberOfComponents(1);
  this->CellPartitionDepths = vtkFloatArray::New();
  this->CellPartitionDepths->SetNumberOfComponents(1);

  this->ToSort = new vtkPCellCenterDepthSortStack;

  this->Timer = vtkTimerLog::New();
}

vtkPCellCenterDepthSort::~vtkPCellCenterDepthSort()
{
  this->SortedCells->Delete();
  this->SortedCellPartition->Delete();
  this->CellCenters->Delete();
  this->CellDepths->Delete();
  this->CellPartitionDepths->Delete();

  delete this->ToSort;
}

void vtkPCellCenterDepthSort::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

float *vtkPCellCenterDepthSort::ComputeProjectionVector()
{
  vtkDebugMacro("ComputeProjectionVector");

  if (this->Camera == NULL)
  {
    vtkErrorMacro("Must set camera before sorting cells.");
    static float v[3] = { 0.0, 0.0, 0.0};
    return v;
  }

  double focalPoint[4];
  double position[4];

  this->Camera->GetFocalPoint(focalPoint);  focalPoint[3] = 1.0;
  this->Camera->GetPosition(position);  position[3] = 1.0;

  this->InverseModelTransform->MultiplyPoint(focalPoint, focalPoint);
  this->InverseModelTransform->MultiplyPoint(position, position);

  static float vector[3];
  if (this->Direction == vtkPVisibilitySort::BACK_TO_FRONT)
  {
    // Sort back to front.
    vector[0] = position[0] - focalPoint[0];
    vector[1] = position[1] - focalPoint[1];
    vector[2] = position[2] - focalPoint[2];
  }
  else
  {
    // Sort front to back.
    vector[0] = focalPoint[0] - position[0];
    vector[1] = focalPoint[1] - position[1];
    vector[2] = focalPoint[2] - position[2];
  }

  vtkDebugMacro("Returning: " << vector[0] << ", " << vector[1] << ", "
                << vector[2]);

  return vector;
}

void vtkPCellCenterDepthSort::ComputeCellCenters()
{
  vtkIdType numcells = this->Input->GetNumberOfCells();
  this->CellCenters->SetNumberOfTuples(numcells);

  float *center = this->CellCenters->GetPointer(0);
  double dcenter[3];
  double *weights = new double[this->Input->GetMaxCellSize()];  //Dummy array.

  for (vtkIdType i = 0; i < numcells; i++)
  {
    vtkCell *cell = this->Input->GetCell(i);
    double pcenter[3];
    int subId;
    subId = cell->GetParametricCenter(pcenter);
    cell->EvaluateLocation(subId, pcenter, dcenter, weights);
    center[0] = dcenter[0]; center[1] = dcenter[1]; center[2] = dcenter[2];
    center += 3;
  }

  delete[] weights;
}

void vtkPCellCenterDepthSort::ComputeDepths()
{
  float *vector = this->ComputeProjectionVector();
  vtkIdType numcells = this->Input->GetNumberOfCells();

  float *center = this->CellCenters->GetPointer(0);
  float *depth = this->CellDepths->GetPointer(0);
  for (vtkIdType i = 0; i < numcells; i++)
  {
    *(depth++) = vtkMath::Dot(center, vector);
    center += 3;
  }
}

void vtkPCellCenterDepthSort::InitTraversal()
{
  vtkDebugMacro("InitTraversal");

  vtkIdType numcells = this->Input->GetNumberOfCells();

  if (   (this->LastSortTime < this->Input->GetMTime())
      || (this->LastSortTime < this->MTime) )
  {
    vtkDebugMacro("Building cell centers array.");

    // Data may have changed.  Recompute cell centers.
    this->ComputeCellCenters();
    this->CellDepths->SetNumberOfTuples(numcells);
    this->SortedCells->SetNumberOfTuples(numcells);
  }

  vtkDebugMacro("Filling SortedCells to initial values.");
  vtkIdType *id = this->SortedCells->GetPointer(0);
  for (vtkIdType i = 0; i < numcells; i++)
  {
    *(id++) = i;
  }

  vtkDebugMacro("Calculating depths.");
  this->ComputeDepths();

  while (!this->ToSort->Stack.empty()) this->ToSort->Stack.pop();
  this->ToSort->Stack.push(vtkIdPair(0, numcells));

  this->LastSortTime.Modified();
}

vtkIdTypeArray *vtkPCellCenterDepthSort::GetNextCells()
{
  if (this->ToSort->Stack.empty())
  {
    // Already sorted and returned everything.
    return NULL;
  }
  vtkIdType *cellIds = this->SortedCells->GetPointer(0);
  float *cellDepths = this->CellDepths->GetPointer(0);
  vtkIdPair partition;

  partition = this->ToSort->Stack.top();  this->ToSort->Stack.pop();
  while (partition.second - partition.first > this->MaxCellsReturned)
  {
    vtkIdType left = partition.first;
    vtkIdType right = partition.second - 1;
    float pivot = cellDepths[static_cast<vtkIdType>(
                               vtkMath::Random(left, right))];
    while (left <= right)
    {
      while ((left <= right) && (cellDepths[left] < pivot)) left++;
      while ((left <= right) && (cellDepths[right] > pivot)) right--;

      if (left > right) break;

      std::swap(cellIds[left], cellIds[right]);
      std::swap(cellDepths[left], cellDepths[right]);

      left++;  right--;
    }

    this->ToSort->Stack.push(vtkIdPair(left, partition.second));
    partition.second = left;
  }

  if (partition.second <= partition.first)
  {
    // Got a partition of zero size.  Just recurse to get the next one.
    return this->GetNextCells();
  }

  vtkIdType firstcell = partition.first;
  vtkIdType numcells = partition.second - partition.first;

  this->SortedCellPartition->SetArray(cellIds + firstcell, numcells, 1);
  this->SortedCellPartition->SetNumberOfTuples(numcells);
  this->CellPartitionDepths->SetArray(cellDepths + firstcell, numcells, 1);
  this->CellPartitionDepths->SetNumberOfTuples(numcells);

  vtkSortDataArray::Sort(this->CellPartitionDepths, this->SortedCellPartition);
  return this->SortedCellPartition;
}


void vtkPCellCenterDepthSort::GetCells()
{
	this->ToSort->Vector.clear();
	this->ToSort->POSVector.clear();
	for (int i = 0; i < this->ToSort->IDVector.size(); ++i)
	{
		this->ToSort->IDVector[i]->Delete();
	}
	this->ToSort->IDVector.clear();

	vtkIdType *cellIds = this->SortedCells->GetPointer(0);
	float *cellDepths = this->CellDepths->GetPointer(0);
	vtkIdPair partition;

	this->Timer->StartTimer();

	while (!this->ToSort->Stack.empty())
	{
		partition = this->ToSort->Stack.top(); this->ToSort->Stack.pop();
		while (partition.second - partition.first > this->MaxCellsReturned)
		{
			vtkIdType left = partition.first;
			vtkIdType right = partition.second - 1;
			float pivot = cellDepths[static_cast<vtkIdType>(
				vtkMath::Random(left, right))];
			while (left <= right)
			{
				while ((left <= right) && (cellDepths[left] < pivot)) left++;
				while ((left <= right) && (cellDepths[right] > pivot)) right--;

				if (left > right) break;

				std::swap(cellIds[left], cellIds[right]);
				std::swap(cellDepths[left], cellDepths[right]);

				left++;  right--;
			}

			this->ToSort->Stack.push(vtkIdPair(left, partition.second));
			partition.second = left;
		}

		if (partition.second > partition.first)
		{
			this->ToSort->Vector.push_back(partition);
		}
	}
	size_t size = this->ToSort->Vector.size();
	this->ToSort->IDVector.resize(size);
	this->ToSort->POSVector.resize(size);

	this->Timer->StopTimer();
	cout << "partical:\t" << this->Timer->GetElapsedTime() << endl;

	this->Timer->StartTimer();

#pragma omp parallel for num_threads(this->Thread)//设置线程数
	for (int i = 0; i < size; i++)
	{
		vtkIdType firstcell = this->ToSort->Vector[i].first;
		vtkIdType numcells = this->ToSort->Vector[i].second - firstcell;
		this->ToSort->IDVector[i]=vtkIdTypeArray::New();
		this->ToSort->IDVector[i]->SetArray(cellIds + firstcell, numcells, 1);
		this->ToSort->IDVector[i]->SetNumberOfTuples(numcells);
		this->ToSort->POSVector[i] = firstcell;
		vtkFloatArray *cellpartitiondepths = vtkFloatArray::New();
		cellpartitiondepths->SetArray(cellDepths + firstcell, numcells, 1);
		cellpartitiondepths->SetNumberOfTuples(numcells);
		vtkSortDataArray::Sort(cellpartitiondepths, this->ToSort->IDVector[i]);
		cellpartitiondepths->Delete();
	}
	this->Timer->StopTimer();
	cout << "global: \t" << this->Timer->GetElapsedTime() << endl;
}

std::vector<vtkIdTypeArray*> vtkPCellCenterDepthSort::GetIDVector()
{
	return this->ToSort->IDVector;
}

std::vector<int> vtkPCellCenterDepthSort::GetPOSVector()
{
	return this->ToSort->POSVector;
}

void vtkPCellCenterDepthSort::SetNumberOfThread(int thread)
{
	this->Thread = thread;
}