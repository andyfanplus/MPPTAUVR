#ifndef vtkPRenderingVolumeOpenGL2ObjectFactory_h
#define vtkPRenderingVolumeOpenGL2ObjectFactory_h

class vtkObjectFactory;
class  vtkPRenderingVolumeOpenGL2ObjectFactory : public vtkObjectFactory
{
public:
	vtkPRenderingVolumeOpenGL2ObjectFactory();
	static vtkPRenderingVolumeOpenGL2ObjectFactory* New()
	{
		vtkPRenderingVolumeOpenGL2ObjectFactory *f = new vtkPRenderingVolumeOpenGL2ObjectFactory;
		f->InitializeObjectBase();
		return f;
	}
	const char* GetVTKSourceVersion() override { return ""; }
	const char* GetDescription() override { return "A fine Test Factory"; }

protected:
	vtkPRenderingVolumeOpenGL2ObjectFactory(const vtkPRenderingVolumeOpenGL2ObjectFactory&) = delete;
	vtkPRenderingVolumeOpenGL2ObjectFactory& operator=(const vtkPRenderingVolumeOpenGL2ObjectFactory&) = delete;
};

#endif // !vtkPRenderingVolumeOpenGL2ObjectFactory_h

