#include <vtkSmartPointer.h>
//#include <vtkHAVSVolumeMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetTriangleFilter.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
#include <vtkVolumeProperty.h>
#include <vtkVolume.h>
#include <vtkCamera.h>
#include <vtkStdString.h>
#include<vtkUnstructuredGrid.h>
#include<vtkUnstructuredGridReader.h>
#include<vtkStructuredGridReader.h>
#include<vtkStructuredGrid.h>
#include<vtkPointData.h>
#include<vtkUnstructuredGridVolumeRayCastMapper.h>
#include<vtkUnstructuredGridVolumeZSweepMapper.h>

#include<vtkMultiBlockDataSet.h>
#include"vtkTecplotReaderBinary.h"

#include<vtkPlot3DMetaReader.h>
#include<vtkMultiBlockPLOT3DReader.h>

#include<vtkObjectFactory.h>
#include"vtkPProjectedTetrahedraMapper.h"
#include"vtkOpenGLPProjectedTetrahedraMapper.h"
#include"vtkPRenderingVolumeOpenGL2ObjectFactory.h"
#include <omp.h>
#include<iostream>
using namespace std;


#define VTK_CREATE(type,name) vtkSmartPointer<type> name=vtkSmartPointer<type>::New()
int main(int argc, char*argv[])
{
	vtkPRenderingVolumeOpenGL2ObjectFactory *f = vtkPRenderingVolumeOpenGL2ObjectFactory::New();
	vtkObjectFactory::RegisterFactory(f);
	f->Delete();

	VTK_CREATE(vtkRenderer, renderer);
	VTK_CREATE(vtkRenderWindow, renwin);

	renwin->AddRenderer(renderer);
	VTK_CREATE(vtkRenderWindowInteractor, iren);
	iren->SetRenderWindow(renwin);


	VTK_CREATE(vtkStructuredGridReader, reader);
	reader->SetFileName("../block71blk.vtk");
	reader->Update();


	VTK_CREATE(vtkDataSetTriangleFilter, filter);
	filter->SetInputData(reader->GetOutput());
	filter->Update();



	double range[2];
	reader->GetOutput()->GetScalarRange(range);
	double mid = (range[0] + range[1]) / 2;


	VTK_CREATE(vtkColorTransferFunction, color);
	color->AddRGBPoint(range[0], 1, 0, 0);
	color->AddRGBPoint(mid, 0, 1, 0);
	color->AddRGBPoint(range[1], 0, 0, 1);
	VTK_CREATE(vtkPiecewiseFunction, opacity);
	opacity->AddPoint(range[0], 0);
	opacity->AddPoint(range[1], 1);



	VTK_CREATE(vtkPProjectedTetrahedraMapper, mapper);
	mapper->SetInputConnection(filter->GetOutputPort());

	mapper->SetNumberOfVST(40);
	mapper->SetNumberOfCDT(40);


	VTK_CREATE(vtkVolumeProperty, property);
	property->SetColor(color);
	property->SetScalarOpacity(opacity);

	vtkVolume *volume = vtkVolume::New();
	volume->SetMapper(mapper);
	volume->SetProperty(property);

	renderer->AddVolume(volume);

	renwin->SetSize(500, 500);
	renwin->Render();
	iren->Initialize();
	iren->Start();


	return EXIT_SUCCESS;
}
