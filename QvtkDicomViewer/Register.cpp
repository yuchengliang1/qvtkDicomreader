#include "Register.h"
#include <QFileDialog>
#include <iostream>

#include<head_all.h>//这是一个不可控点

#include "itkGDCMImageIO.h"  
#include "itkImageRegistrationMethodv4.h"
#include "itkImageRegistrationMethod.h"
#include "itkTranslationTransform.h"
#include "itkMeanSquaresImageToImageMetricv4.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkRegularStepGradientDescentOptimizerv4.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRegularStepGradientDescentBaseOptimizer.h"
#include <vtkObjectFactory.h>
#include <vtkLookupTable.h>
#include <vtkImageHistogramStatistics.h>
#include <QMessageBox>

void wheelCancelInteractorStyle::OnMouseWheelForward()
{
}
void wheelCancelInteractorStyle::OnMouseWheelBackward()
{
}
vtkStandardNewMacro(wheelCancelInteractorStyle);

class vtkImageInteractionCallback : public vtkCommand
{
public:
	static vtkImageInteractionCallback* New() { return new vtkImageInteractionCallback; }

	vtkImageInteractionCallback()
	{
		this->Viewer = nullptr;
	}

	void SetImageViewer(myVtkViewer* viewer) { this->Viewer = viewer; }
	myVtkViewer* GetImageViewer() { return this->Viewer; }

	void Execute(vtkObject*, unsigned long event, void*) override
	{
		if (event == vtkCommand::MouseWheelForwardEvent)
		{
			myVtkViewer* viewer = this->GetImageViewer();
			auto sliceIndex = viewer->GetSlice();
			viewer->SetSlice(sliceIndex + 1);
			viewer->Render();
		}
		else if (event == vtkCommand::MouseWheelBackwardEvent)
		{
			myVtkViewer* viewer = this->GetImageViewer();
			auto sliceIndex = viewer->GetSlice();
			viewer->SetSlice(sliceIndex - 1);
			viewer->Render();
		}
	}

private:
	myVtkViewer* Viewer;
};

/*
 * 构造
 */
Register::Register(QWidget *parent)
	: QWidget(parent)
{
	ui.setupUi(this);
	for(int i=0;i<4;i++ )
	{
		actor[i] = vtkSmartPointer<vtkImageActor>::New();
		renderer[i] = vtkSmartPointer<vtkRenderer>::New();
		renderWindowInteractor[i] = vtkSmartPointer<vtkRenderWindowInteractor>::New();
		style[i] = vtkSmartPointer<wheelCancelInteractorStyle>::New();
		viewer[i] = vtkSmartPointer<myVtkViewer>::New();
		connector[i] = ConnectorType::New();
		// colors[i] = vtkSmartPointer<vtkImageMapToColors>::New();
	}
	m_output_widgets[0]=ui.qvtkWidget_Registration_UL;
	m_output_widgets[1]=ui.qvtkWidget_Registration_UR;
	m_output_widgets[2]=ui.qvtkWidget_Registration_DL;
	m_output_widgets[3]=ui.qvtkWidget_Registration_DR;
	m_CurrentRegFunc = RegFunc_Translation;//默认
	fixedImageReader = FixedImageReaderType::New();
	movingImageReader = MovingImageReaderType::New();
	defaultImageReader= FixedImageReaderType::New();
	/*
	 * 加载一张默认图片刷掉白屏
	 */

	//defaultImageReader->SetFileName("./Resources/default.png");
	//for(int i=0;i<4;i++)
	//{
	//	connector[i]->SetInput(defaultImageReader->GetOutput());
	//	connector[i]->Update();
	//	actor[i]->GetMapper()->SetInputData(connector[i]->GetOutput());
	//	renderer[i]->AddActor(actor[i]);
	//	m_output_widgets[i]->GetRenderWindow()->AddRenderer(renderer[i]);
	//	renderWindowInteractor[i]->SetRenderWindow(m_output_widgets[i]->GetRenderWindow());
	//	renderWindowInteractor[i]->SetInteractorStyle(style[i]);
	//	renderWindowInteractor[i]->Initialize();
	//	m_output_widgets[i]->GetRenderWindow()->Render();
	//}
	//renderWindowInteractor[0]->Start();
	//renderWindowInteractor[1]->Start();
	//renderWindowInteractor[2]->Start();
	//renderWindowInteractor[3]->Start();
}
/*
 * 析构
 */
Register::~Register()
{
}

/*
	3D仿射变换
*/
void Register::AffineTransformReg(FixedImageReaderType::Pointer _fixedImageReader, MovingImageReaderType::Pointer _movingImageReader, double steplength, unsigned int maxNumberOfIterations)
{
	constexpr unsigned int Dimension = 3;
	using PixelType = float;

	typedef itk::Image<PixelType, Dimension> FixedImageType;
	typedef itk::Image<PixelType, Dimension> MovingImageType;
	typedef itk::AffineTransform<double, Dimension> TransformType;
	typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
	typedef itk::MeanSquaresImageToImageMetric<FixedImageType, MovingImageType> MetricType;
	typedef itk::LinearInterpolateImageFunction<MovingImageType, double> InterpolatorType;
	typedef itk::ImageRegistrationMethod<FixedImageType, MovingImageType> RegistrationType;

	MetricType::Pointer metric = MetricType::New();
	OptimizerType::Pointer optimizer = OptimizerType::New();
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	RegistrationType::Pointer registration = RegistrationType::New();
	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);
	registration->SetInterpolator(interpolator);

	TransformType::Pointer transform = TransformType::New();
	registration->SetTransform(transform);
	registration->SetFixedImage(_fixedImageReader->GetOutput());
	registration->SetMovingImage(_movingImageReader->GetOutput());
	_fixedImageReader->Update();
	registration->SetFixedImageRegion(_fixedImageReader->GetOutput()->GetBufferedRegion());

	typedef itk::CenteredTransformInitializer<TransformType, FixedImageType, MovingImageType> TransformInitializerType;
	TransformInitializerType::Pointer initializer = TransformInitializerType::New();

	initializer->SetTransform(transform);
	initializer->SetFixedImage(_fixedImageReader->GetOutput());
	initializer->SetMovingImage(_movingImageReader->GetOutput());
	initializer->MomentsOn();
	initializer->InitializeTransform();

	registration->SetInitialTransformParameters(transform->GetParameters());

	double translationScale = 1.0 / 1000.0;

	using OptimizerScalesType = OptimizerType::ScalesType;
	OptimizerScalesType optimizerScales(transform->GetNumberOfParameters());
	optimizerScales[0] = 1.0;
	optimizerScales[1] = 1.0;
	optimizerScales[2] = 1.0;
	optimizerScales[3] = 1.0;
	optimizerScales[4] = 1.0;
	optimizerScales[5] = 1.0;
	optimizerScales[6] = 1.0;
	optimizerScales[7] = 1.0;
	optimizerScales[8] = 1.0;
	optimizerScales[9] = translationScale;
	optimizerScales[10] = translationScale;
	optimizerScales[11] = translationScale;
	optimizer->SetScales(optimizerScales);


	steplength = 0.1;
	maxNumberOfIterations = 2;

	optimizer->SetMaximumStepLength(steplength);
	optimizer->SetMinimumStepLength(0.0001);
	optimizer->SetNumberOfIterations(maxNumberOfIterations);
	optimizer->MinimizeOn();

	try
	{
		registration->Update();
		std::cout << "Optimizer stop condition: "
			<< registration->GetOptimizer()->GetStopConditionDescription()
			<< std::endl;
	}
	catch (const itk::ExceptionObject& err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		return;
	}

	OptimizerType::ParametersType finalParameters = registration->GetLastTransformParameters();
	const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
	const double bestValue = optimizer->GetValue();

	using ResampleFilterType = itk::ResampleImageFilter<MovingImageType, FixedImageType>;

	TransformType::Pointer finalTransform = TransformType::New();
	finalTransform->SetParameters(finalParameters);
	finalTransform->SetFixedParameters(transform->GetFixedParameters());

	auto resampler = ResampleFilterType::New();
	resampler->SetTransform(finalTransform);
	resampler->SetInput(_movingImageReader->GetOutput());

	FixedImageType::Pointer fixedImage = _fixedImageReader->GetOutput();
	resampler->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
	resampler->SetOutputOrigin(fixedImage->GetOrigin());
	resampler->SetOutputSpacing(fixedImage->GetSpacing());
	resampler->SetOutputDirection(fixedImage->GetDirection());
	resampler->SetDefaultPixelValue(0);

	using OutputPixelType = float;
	using OutputImageType = itk::Image<OutputPixelType, Dimension>;
	using CastFilterType = itk::CastImageFilter<FixedImageType, OutputImageType>;
	using WriterType = itk::ImageFileWriter<OutputImageType>;

	auto caster = CastFilterType::New();
	caster->SetInput(resampler->GetOutput());

	using RescalerType = itk::RescaleIntensityImageFilter<FixedImageType, OutputImageType>;
	auto intensityRescaler = RescalerType::New();
	intensityRescaler->SetInput(caster->GetOutput());
	intensityRescaler->SetOutputMinimum(0);
	intensityRescaler->SetOutputMaximum(255);
	connector[2]->SetInput(intensityRescaler->GetOutput());
	//connector[2]->SetInput(caster->GetOutput());

	auto writer = WriterType::New();
	writer->SetFileName("registration.nii.gz");
	writer->SetInput(intensityRescaler->GetOutput());
	writer->Update();

	using DifferenceFilterType = itk::SubtractImageFilter<FixedImageType, FixedImageType, FixedImageType>;
	auto difference = DifferenceFilterType::New();
	difference->SetInput1(_fixedImageReader->GetOutput());
	difference->SetInput2(resampler->GetOutput());
	auto intensityRescaler2 = RescalerType::New();
	intensityRescaler2->SetInput(difference->GetOutput());
	intensityRescaler2->SetOutputMinimum(0);
	intensityRescaler2->SetOutputMaximum(255);
	connector[3]->SetInput(intensityRescaler2->GetOutput());

	auto writer2 = WriterType::New();
	writer2->SetFileName("checker.nii.gz");
	writer2->SetInput(intensityRescaler2->GetOutput());
	writer2->Update();
	
	setWindowTitle(QStringLiteral("配准工具:计算完成!"));

	//vtkSmartPointer<vtkLookupTable> greyLut = vtkSmartPointer<vtkLookupTable>::New();
	//greyLut->SetNumberOfTableValues(256);
	//greyLut->SetTableRange(0, 255);
	//greyLut->SetHueRange(0.0, 0.0);
	//greyLut->SetSaturationRange(0.0, 0.0);
	//greyLut->SetValueRange(0.0, 1.0);
	//greyLut->SetAlphaRange(1.0, 1.0);
	//greyLut->Build();

	for (int i = 2; i < 4; i++)
	{
		connector[i]->Update();

		//colors[i]->SetInputData(connector[i]->GetOutput());
		//colors[i]->SetLookupTable(greyLut);
		//colors[i]->Update();

		viewer[i]->SetInputData(connector[i]->GetOutput());
		viewer[i]->SetSize(340, 255);

		vtkSmartPointer<vtkImageHistogramStatistics> stats = vtkSmartPointer<vtkImageHistogramStatistics>::New();
		stats->SetAutoRangePercentiles(0.1, 99.9);
		stats->SetAutoRangeExpansionFactors(0.0, 0.0);
		stats->SetInputData(connector[i]->GetOutput());
		stats->Update();
		double* intensityRange = stats->GetAutoRange();
		double window = intensityRange[1] - intensityRange[0];
		double level = 0.5 * (intensityRange[0] + intensityRange[1]);

		viewer[i]->SetColorLevel(level); // 窗位
		viewer[i]->SetColorWindow(window); // 窗宽
		viewer[i]->SetSliceOrientationToXY();
		viewer[i]->SetSlice(0);

		viewer[i]->SetRenderer(renderer[i]);
		// m_output_widgets[0]->GetRenderWindow()->AddRenderer(renderer[0]);
		viewer[i]->SetRenderWindow(m_output_widgets[i]->GetRenderWindow());
		viewer[i]->Render();

		// renderWindowInteractor[0]->SetRenderWindow(m_output_widgets[0]->GetRenderWindow());
		vtkSmartPointer<vtkImageInteractionCallback> callback = vtkSmartPointer<vtkImageInteractionCallback>::New();
		callback->SetImageViewer(viewer[i]);

		renderWindowInteractor[i]->AddObserver(vtkCommand::MouseWheelForwardEvent, callback);
		renderWindowInteractor[i]->AddObserver(vtkCommand::MouseWheelBackwardEvent, callback);

		renderWindowInteractor[i]->SetInteractorStyle(style[i]);
		viewer[i]->SetupInteractor(renderWindowInteractor[i]);
	}
	for (int i = 2; i < 4; i++)
	{
		m_output_widgets[i]->GetRenderWindow()->Render();
	}
	//下边的这个有阻塞,不要放在循环里
	renderWindowInteractor[2]->Start();
	renderWindowInteractor[3]->Start();
}

/*
 * 刷新图片
 */
void Register::updateOutputImage()
{
	setWindowTitle(QStringLiteral("配准工具:计算完成!"));

	vtkSmartPointer<vtkLookupTable> greyLut = vtkSmartPointer<vtkLookupTable>::New();
	greyLut->SetNumberOfTableValues(256);
	greyLut->SetTableRange(0, 255);
	greyLut->SetHueRange(0.0, 0.0);
	greyLut->SetSaturationRange(0.0, 0.0);
	greyLut->SetValueRange(0.0, 1.0);
	greyLut->SetAlphaRange(1.0, 1.0);
	greyLut->Build();

	for (int i = 2; i<4; i++)
	{
		connector[i]->Update();

		//colors[i]->SetInputData(connector[i]->GetOutput());
		//colors[i]->SetLookupTable(greyLut);
		//colors[i]->Update();

		//viewer[i]->SetInputData(colors[i]->GetOutput());
		viewer[i]->SetSize(340, 255);

		vtkSmartPointer<vtkImageHistogramStatistics> stats = vtkSmartPointer<vtkImageHistogramStatistics>::New();
		stats->SetAutoRangePercentiles(0.1, 99.9);
		stats->SetAutoRangeExpansionFactors(0.0, 0.0);
		stats->SetInputData(connector[i]->GetOutput());
		stats->Update();
		double* intensityRange = stats->GetAutoRange();
		double window = intensityRange[1] - intensityRange[0];
		double level = 0.5 * (intensityRange[0] + intensityRange[1]);

		viewer[i]->SetColorLevel(62); // 窗位
		viewer[i]->SetColorWindow(124); // 窗宽
		viewer[i]->SetSliceOrientationToXY();
		viewer[i]->SetSlice(0);

		viewer[i]->SetRenderer(renderer[i]);
		// m_output_widgets[0]->GetRenderWindow()->AddRenderer(renderer[0]);
		viewer[i]->SetRenderWindow(m_output_widgets[i]->GetRenderWindow());
		viewer[i]->Render();

		// renderWindowInteractor[0]->SetRenderWindow(m_output_widgets[0]->GetRenderWindow());
		vtkSmartPointer<vtkImageInteractionCallback> callback = vtkSmartPointer<vtkImageInteractionCallback>::New();
		callback->SetImageViewer(viewer[i]);

		renderWindowInteractor[i]->AddObserver(vtkCommand::MouseWheelForwardEvent, callback);
		renderWindowInteractor[i]->AddObserver(vtkCommand::MouseWheelBackwardEvent, callback);

		renderWindowInteractor[i]->SetInteractorStyle(style[i]);
		viewer[i]->SetupInteractor(renderWindowInteractor[i]);
	}
	for (int i = 2; i<4; i++)
	{
		m_output_widgets[i]->GetRenderWindow()->Render();
	}
	//下边的这个有阻塞,不要放在循环里
	renderWindowInteractor[2]->Start();
	renderWindowInteractor[3]->Start();
}

/*
 * 选择基准图片
 */
void Register::OnSelectImageFix()
{
	QString fileName = QFileDialog::getOpenFileName(this, QStringLiteral("选择基准图像"), NULL, tr("*.*"));
	ui.lineEdit_ImageFix->setText(fileName);
	if (fileName.isEmpty())
	{
		QMessageBox::information(this, QStringLiteral("提示"), QStringLiteral("请选择文件图片!"));
	}
	else
	{
		try
		{
			fixedImageReader->SetFileName(fileName.toStdString());
			fixedImageReader->Update();

			connector[0]->SetInput(fixedImageReader->GetOutput());//显示输入1
			connector[0]->Update();

			viewer[0]->SetInputData(connector[0]->GetOutput());

			vtkSmartPointer<vtkImageHistogramStatistics> stats0 = vtkSmartPointer<vtkImageHistogramStatistics>::New();
			stats0->SetAutoRangePercentiles(0.1, 99.9);
			stats0->SetAutoRangeExpansionFactors(0.0, 0.0);
			stats0->SetInputData(connector[0]->GetOutput());
			stats0->Update();
			double* intensityRange0 = stats0->GetAutoRange();
			double window = intensityRange0[1] - intensityRange0[0];
			double level = 0.5 * (intensityRange0[0] + intensityRange0[1]);

			viewer[0]->SetColorLevel(level); // 窗位
			viewer[0]->SetColorWindow(window); // 窗宽

			viewer[0]->SetSize(340, 255);
			viewer[0]->SetSliceOrientationToXY();
			viewer[0]->SetSlice(0);
			viewer[0]->SetRenderer(renderer[0]);
			// m_output_widgets[0]->GetRenderWindow()->AddRenderer(renderer[0]);
			viewer[0]->SetRenderWindow(m_output_widgets[0]->GetRenderWindow());
			viewer[0]->Render();

			// renderWindowInteractor[0]->SetRenderWindow(m_output_widgets[0]->GetRenderWindow());
			vtkSmartPointer<vtkImageInteractionCallback> callback = vtkSmartPointer<vtkImageInteractionCallback>::New();
			callback->SetImageViewer(viewer[0]);

			renderWindowInteractor[0]->AddObserver(vtkCommand::MouseWheelForwardEvent, callback);
			renderWindowInteractor[0]->AddObserver(vtkCommand::MouseWheelBackwardEvent, callback);

			renderWindowInteractor[0]->SetInteractorStyle(style[0]);
			viewer[0]->SetupInteractor(renderWindowInteractor[0]);
			renderWindowInteractor[0]->Start();
		}
		catch (itk::ExceptionObject & err)
		{
			QString errorMsg= err.GetDescription();
			QMessageBox::information(this, QStringLiteral("发生严重错误!"), errorMsg);
		}
		
	}
}
/*
 * 选择待配准图片
 */
void Register::OnSelectImageMove()
{
	QString fileName = QFileDialog::getOpenFileName(this, QStringLiteral("选择待配准图像"), NULL, tr("*.*"));
	ui.lineEdit_ImageMove->setText(fileName);
	if (fileName.isEmpty())
	{
		QMessageBox::information(this, QStringLiteral("提示"), QStringLiteral("请选择文件图片!"));
	}
	else
	{
		try
		{
			movingImageReader->SetFileName(fileName.toStdString());
			movingImageReader->Update();
			connector[1]->SetInput(movingImageReader->GetOutput());//显示输入1
			connector[1]->Update();

			viewer[1]->SetInputData(connector[1]->GetOutput());

			viewer[1]->SetSize(340, 255);
			viewer[1]->SetColorLevel(770); // 窗位
			viewer[1]->SetColorWindow(1541); // 窗宽
			viewer[1]->SetSliceOrientationToXY();
			viewer[1]->SetSlice(0);

			viewer[1]->SetRenderer(renderer[1]);
			// m_output_widgets[0]->GetRenderWindow()->AddRenderer(renderer[0]);
			viewer[1]->SetRenderWindow(m_output_widgets[1]->GetRenderWindow());
			viewer[1]->Render();

			// renderWindowInteractor[0]->SetRenderWindow(m_output_widgets[0]->GetRenderWindow());
			vtkSmartPointer<vtkImageInteractionCallback> callback = vtkSmartPointer<vtkImageInteractionCallback>::New();
			callback->SetImageViewer(viewer[1]);

			renderWindowInteractor[1]->AddObserver(vtkCommand::MouseWheelForwardEvent, callback);
			renderWindowInteractor[1]->AddObserver(vtkCommand::MouseWheelBackwardEvent, callback);

			renderWindowInteractor[1]->SetInteractorStyle(style[1]);
			viewer[1]->SetupInteractor(renderWindowInteractor[1]);
			renderWindowInteractor[1]->Start();
		}
		catch (itk::ExceptionObject & err)
		{
			QString errorMsg = err.GetDescription();
			QMessageBox::information(this, QStringLiteral("发生严重错误!"), errorMsg);
		}

	}
}
/*
 * 开始计算
 */
void Register::OnButtonOk()
{
	setWindowTitle(QStringLiteral("配准工具:正在计算..."));
	switch (m_CurrentRegFunc)
	{
	case RegFunc_Translation:
		AffineTransformReg(fixedImageReader, movingImageReader);
		break;
	case RegFunc_CenteredSimilarity:
		AffineTransformReg(fixedImageReader, movingImageReader);
		break;
	case RegFunc_Affine:
		AffineTransformReg(fixedImageReader, movingImageReader);
		break;
	case RegFunc_Multi:
		AffineTransformReg(fixedImageReader, movingImageReader);
		break;
	default:
		break;
	}
	//this->close();
}
/*
 * 退出
 */
void Register::OnButtonCancel()
{
	this->close();
}
/*
 * 选择平移变换
 */
void Register::OnSelectTranslation()
{
	m_CurrentRegFunc = RegFunc_Translation;
}
/*
 * 选择中心对称二维变换
 */
void Register::OnSelectCenteredSimilarity()
{
	m_CurrentRegFunc = RegFunc_CenteredSimilarity;
}
/*
 * 选择仿射变换
 */
void Register::OnSelectAffine()
{
	m_CurrentRegFunc = RegFunc_Affine;
}
/*
 * 选择膜法变换
 */
void Register::OnSelectMulit()
{
	m_CurrentRegFunc = RegFunc_Multi;
}
