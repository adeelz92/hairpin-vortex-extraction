// In this version, we use both geometric and physical criteria to simplify regions.
// In this version, we try to simplify the region using vtkExtractSurface filter.
// In this version, we use vtkPoissonReconstruction to extract surfaces.
// In this version, we extract the skeleton using the compiled exe file of the GCAL MCF Skeletonization.
// In this version, we smooth-out the skeleton after extraction.
// In this version, we assign a unique region id to the skeleton after extraction.
// This version saves the surface data alongwith skeleton.

// Skeleton_Extraction.exe "Path\To\Dataset_SimpleRegions.vtk" "Path\To\Dataset_Skeleton.vtk" "Path\To\Dataset_Surface.vtk"
#include <vtkActor.h>
#include <vtkAppendPolyData.h>
#include <vtkAppendFilter.h>
#include <vtkCamera.h>
#include <vtkCleanPolyData.h>
#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>
#include <vtkExtractSurface.h>
#include <vtkGaussianKernel.h>
#include <vtkGeometryFilter.h>
#include <vtkFillHolesFilter.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkjsoncpp/json/json.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPCANormalEstimation.h>
#include <vtkPointData.h>
#include <vtkPointDataToCellData.h>
#include <vtkPointInterpolator.h>
#include <vtkPoissonReconstruction.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataPointSampler.h>
#include <vtkPolyLine.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkSMPTools.h>
#include <vtkSTLWriter.h>
#include <vtkThreshold.h>
#include <vtkUnstructuredGrid.h>

#include <stdlib.h>
#include <chrono>
#include <conio.h>
#include <iostream>
#include <fstream>
#include <random>

template <typename T>
struct Region
{
	int RegionId;
	vtkSmartPointer<T> RegionData;

	Region(int regionId, T* regionData)
	{
		this->RegionId = regionId;
		this->RegionData = vtkSmartPointer<T>::New();
		this->RegionData->DeepCopy(regionData);
	}

	void updateData(T* regionData)
	{
		this->RegionData->DeepCopy(regionData);
	}
};

void extractSurface(vtkDataSet* dataset, vtkPolyData* output)
{
	/**/
	if (dataset->GetNumberOfPoints() == 0 || dataset->GetNumberOfCells() == 0)
	{
		output = nullptr;
		return;
	}
	vtkSmartPointer<vtkGeometryFilter> surfaceFilter = vtkSmartPointer<vtkGeometryFilter>::New();
	surfaceFilter->SetInputData(dataset);
	surfaceFilter->Update();

	if (surfaceFilter->GetOutput()->GetNumberOfPoints() == 0 || surfaceFilter->GetOutput()->GetNumberOfCells() == 0)
	{
		output = nullptr;
		return;
	}

	vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smoothFilter->SetInputConnection(surfaceFilter->GetOutputPort());
	smoothFilter->SetNumberOfIterations(10);
	smoothFilter->SetRelaxationFactor(0.25);
	smoothFilter->FeatureEdgeSmoothingOff();
	smoothFilter->BoundarySmoothingOn();
	smoothFilter->Update();

	if (smoothFilter->GetOutput()->GetNumberOfPoints() == 0 || smoothFilter->GetOutput()->GetNumberOfCells() == 0)
	{
		output = nullptr;
		return;
	}

	/*
	int cellsBefore = smoothFilter->GetOutput()->GetNumberOfCells();
	vtkSmartPointer<vtkFillHolesFilter> fillHoles = vtkSmartPointer<vtkFillHolesFilter>::New();
	fillHoles->SetInputConnection(smoothFilter->GetOutputPort());
	fillHoles->SetHoleSize(smoothFilter->GetOutput()->GetLength()/4);
	fillHoles->Update();

	// Make the triangle winding order consistent
	vtkNew<vtkPolyDataNormals> normals;
	normals->SetInputData(fillHoles->GetOutput());
	normals->ConsistencyOn();
	normals->SplittingOff();
	normals->Update();
	int cellsAfter = normals->GetOutput()->GetNumberOfCells();
	if (cellsBefore != cellsAfter)
	{
		cout << "Number of cells before: " << cellsBefore << " After: " << cellsAfter << endl;
	}
	*/

	vtkPolyData* surfaceData = smoothFilter->GetOutput();
	output->DeepCopy(surfaceData);
	return;

	double bounds[6];
	surfaceData->GetBounds(bounds);
	double range[3];
	for (int i = 0; i < 3; ++i)
	{
		range[i] = bounds[2 * i + 1] - bounds[2 * i];
	}
	std::cout << "# of original points: " << surfaceData->GetNumberOfPoints() << std::endl;
	vtkNew<vtkPolyDataPointSampler> sample;
	sample->SetInputData(surfaceData);
	// std::min(std::min(range[0], range[1]), range[2]) / 5
	sample->SetDistance(1.5);
	sample->Update();
	std::cout << "# of points after sampling: " << sample->GetOutput()->GetNumberOfPoints() << std::endl;
	surfaceData = sample->GetOutput();

	vtkSmartPointer<vtkPoissonReconstruction> surface = vtkSmartPointer<vtkPoissonReconstruction>::New();
	surface->SetDepth(7);
	//surface->SetKernelDepth(7);

	int sampleSize = surfaceData->GetNumberOfPoints() * 0.0005;
	if (sampleSize < 10)
	{
		sampleSize = 10;
	}
	std::cout << "Sample size is: " << sampleSize << std::endl;

	if (surfaceData->GetPointData()->GetNormals())
	{
		surface->SetInputData(surfaceData);
	}
	else
	{
		vtkSmartPointer<vtkPCANormalEstimation> normals = vtkSmartPointer<vtkPCANormalEstimation>::New();
		normals->SetInputData(surfaceData);
		normals->SetSampleSize(sampleSize);
		normals->SetNormalOrientationToGraphTraversal();
		normals->FlipNormalsOn();
		surface->SetInputConnection(normals->GetOutputPort());
	}

	/*
	double bounds[6];
	surfaceData->GetBounds(bounds);
	double range[3];
	for (int i = 0; i < 3; ++i)
	{
		range[i] = bounds[2 * i + 1] - bounds[2 * i];
	}

	std::cout << "# of original points: " << surfaceData->GetNumberOfPoints() << std::endl;
	vtkNew<vtkPolyDataPointSampler> sample;
	sample->SetInputData(surfaceData);
	sample->SetDistance(range[0] / 20);
	sample->Update();
	std::cout << "# of points after sampling: " << sample->GetOutput()->GetNumberOfPoints() << std::endl;
	surfaceData = sample->GetOutput();

	surfaceData->GetBounds(bounds);
	for (int i = 0; i < 3; ++i)
	{
		range[i] = bounds[2 * i + 1] - bounds[2 * i];
	}

	std::cout << "Range: " << range[0] << ", " << range[1] << ", " << range[2] << std::endl;
	int dimension = 256;
	double radius;
	radius = std::max(std::max(range[0], range[1]), range[2]) / static_cast<double>(dimension) * 25;
	distance->SetRadius(radius);
	distance->SetDimensions(dimension, dimension, dimension);
	distance->SetBounds(bounds[0] - range[0] * .25, bounds[1] + range[0] * .25,
		bounds[2] - range[1] * .25, bounds[3] + range[1] * .25,
		bounds[4] - range[2] * .25, bounds[5] + range[2] * .25);

	vtkSmartPointer<vtkExtractSurface> surface = vtkSmartPointer<vtkExtractSurface>::New();
	surface->SetInputConnection(distance->GetOutputPort());
	surface->SetRadius(radius * .99);
	// surface->Update();
	*/
	vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
	cleanFilter->SetInputConnection(surface->GetOutputPort());
	cleanFilter->Update();

	output->DeepCopy(cleanFilter->GetOutput());

}

double round_up(double value, int decimal_places)
{
	const double multiplier = std::pow(10.0, decimal_places);
	return std::ceil(value * multiplier) / multiplier;
}

void readSkelInfo(std::string& filename, vtkPolyData* skeleton)
{
	fstream fin;
	fin.open(filename.c_str(), std::ios::in);
	if (!fin.is_open()) throw std::runtime_error("Could not open file");
	std::vector<std::array<double, 3>> values;
	vtkNew<vtkPoints> points;
	vtkNew<vtkCellArray> cells;
	std::string line;
	double val, numPts;
	vtkSmartPointer<vtkPolyLine> polyLine = nullptr;

	int ptIdx = -1;
	int linePtIdx;
	while (std::getline(fin, line))
	{
		std::stringstream ss(line);
		std::vector<double> results;

		while (ss >> val) {
			results.push_back(val);
		}

		if (results.size() == 1)
		{
			linePtIdx = 0;
			numPts = results.at(0);
			// cout << numPts << endl;
			if (polyLine != nullptr)
			{
				cells->InsertNextCell(polyLine);
			}
			polyLine = vtkSmartPointer<vtkPolyLine>::New();
			polyLine->GetPointIds()->SetNumberOfIds(numPts);
		}

		if (results.size() > 1)
		{
			double pt[3] = { round_up(results.at(0),4), round_up(results.at(1),4), round_up(results.at(2),4) };
			/**/
			// cout << results.at(0) << " " << results.at(1) << " " << results.at(2);
			bool found = false;
			int i;
			for (i = 0; i < points->GetNumberOfPoints(); i++)
			{
				double x[3];
				points->GetPoint(i, x);
				x[0] = round_up(x[0], 4); x[1] = round_up(x[1], 4); x[2] = round_up(x[2], 4);
				if (pt[0] == x[0] && pt[1] == x[1] && pt[2] == x[2])
				{
					found = true;
					break;
				}
			}
			// cout << " " << polyLine->GetNumberOfPoints() << endl;
			if (!found)
			{
				ptIdx++;
				points->InsertNextPoint(pt);
				polyLine->GetPointIds()->SetId(linePtIdx, ptIdx);
				linePtIdx++;
			}
			else
			{
				polyLine->GetPointIds()->SetId(linePtIdx, i);
				linePtIdx++;
			}
		}
	}
	cells->InsertNextCell(polyLine);

	// Add the points to the dataset
	skeleton->SetPoints(points);
	// Add the lines to the dataset
	skeleton->SetLines(cells);
}

void extractSkeleton(vtkDataSet* dataset, vtkPolyData* output,
	vtkPolyData* surface_output, int regionId)
{
	vtkNew<vtkNamedColors> colors;

	vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();
	extractSurface(dataset, surface);
	if (surface == nullptr)
	{
		output = nullptr;
		surface_output = nullptr;
		return;
	}
	surface_output->DeepCopy(surface);

	std::string filename = "mesh";
	std::string ext = ".stl";
	std::string mesh_filename_wExt = filename + ext;
	// First write the extracted surface to an STL file
	std::string mesh_filepath = mesh_filename_wExt;

	vtkNew<vtkSTLWriter> stlWriter;
	stlWriter->SetFileName(mesh_filepath.c_str());
	stlWriter->SetInputData(surface);
	stlWriter->SetFileTypeToBinary();
	stlWriter->Write();

	std::string skel_exe_path = "Extract_Skeleton.exe";
	std::string skel_output_path = filename + ".txt";
	std::string command = skel_exe_path + " \"" + mesh_filepath + "\" \"" + skel_output_path + "\"";
	system(command.c_str());

	vtkSmartPointer<vtkPolyData> skeleton = vtkSmartPointer<vtkPolyData>::New();
	readSkelInfo(skel_output_path, skeleton);

	vtkSmartPointer<vtkGaussianKernel> kernel = vtkSmartPointer<vtkGaussianKernel>::New();
	kernel->SetKernelFootprintToNClosest();
	kernel->SetNumberOfPoints(8);

	vtkSmartPointer<vtkPointInterpolator> interpolator = vtkSmartPointer<vtkPointInterpolator>::New();
	interpolator->SetInputData(skeleton);
	interpolator->SetSourceData(dataset);
	interpolator->SetKernel(kernel);
	interpolator->SetNullPointsStrategyToClosestPoint();
	interpolator->Update();
	skeleton = interpolator->GetPolyDataOutput();

	// Add regionIds array
	vtkNew<vtkIdTypeArray> skelRegionIdsArray;
	skelRegionIdsArray->SetName("skeletonRegionIds");
	skelRegionIdsArray->SetNumberOfTuples(skeleton->GetNumberOfCells());
	for (int i = 0; i < skeleton->GetNumberOfCells(); i++)
	{
		skelRegionIdsArray->SetTuple1(i, regionId);
	}
	skeleton->GetCellData()->AddArray(skelRegionIdsArray);

	vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smoothFilter->SetInputData(skeleton);
	smoothFilter->SetNumberOfIterations(10);
	smoothFilter->SetRelaxationFactor(0.5);
	smoothFilter->FeatureEdgeSmoothingOff();
	smoothFilter->BoundarySmoothingOn();
	smoothFilter->Update();
	vtkPolyData* smoothSkeleton = smoothFilter->GetOutput();
	output->DeepCopy(smoothSkeleton);
	/*
	cout << *smoothSkeleton << endl;

	vtkNew<vtkRenderer> ren1;
	ren1->SetBackground(colors->GetColor3d("White").GetData());
	ren1->SetViewport(0, 0, 1, 0.5);
	vtkCamera* camera = ren1->GetActiveCamera();
	vtkNew<vtkRenderer> ren2;
	ren2->SetBackground(colors->GetColor3d("White").GetData());
	ren2->SetViewport(0, 0.5, 1, 1);
	ren2->SetActiveCamera(camera);

	vtkNew<vtkDataSetMapper> skelMapper;
	skelMapper->SetInputData(skeleton);
	skelMapper->ScalarVisibilityOff();
	skelMapper->Update();
	vtkNew<vtkActor> skelActor;
	skelActor->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());
	skelActor->GetProperty()->SetLineWidth(2.0);
	skelActor->GetProperty()->SetOpacity(1);
	skelActor->SetMapper(skelMapper);
	ren2->AddActor(skelActor);

	vtkNew<vtkDataSetMapper> regMapper;
	regMapper->SetInputData(smoothSkeleton);
	regMapper->ScalarVisibilityOff();
	regMapper->Update();
	vtkNew<vtkActor> regActor;
	regActor->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());
	regActor->SetMapper(regMapper);
	regActor->GetProperty()->SetLineWidth(2.0);
	regActor->GetProperty()->SetOpacity(1);
	ren1->AddActor(regActor);

	vtkNew<vtkDataSetMapper> smallRegMapper;
	smallRegMapper->SetInputData(surface);
	smallRegMapper->ScalarVisibilityOff();
	smallRegMapper->Update();
	vtkNew<vtkActor> smallRegActor;
	smallRegActor->GetProperty()->SetColor(colors->GetColor3d("Green").GetData());
	smallRegActor->GetProperty()->SetLineWidth(2.0);
	smallRegActor->GetProperty()->SetOpacity(0.2);
	smallRegActor->SetMapper(smallRegMapper);
	ren2->AddActor(smallRegActor);

	vtkNew<vtkRenderWindow> ren_win;
	ren_win->SetSize(1200, 600);
	ren_win->SetWindowName("Vortex Lines");
	ren_win->AddRenderer(ren1);
	ren_win->AddRenderer(ren2);
	ren1->ResetCamera();

	vtkNew<vtkRenderWindowInteractor> iren;
	vtkNew<vtkInteractorStyleTrackballCamera> irenStyle;
	iren->SetInteractorStyle(irenStyle);
	iren->SetRenderWindow(ren_win);

	ren_win->Render();
	iren->Start();
	*/
}

int main(int argc, char* argv[])
{
	vtkSMPTools::SetBackend("STDThread");
	// vtkSMPTools::Initialize(8);
	int num_threads = vtkSMPTools::GetEstimatedNumberOfThreads();
	cout << "Total threads: " << num_threads << endl;

	vtkNew<vtkNamedColors> colors;

	if (argc < 4)
	{
		std::cerr << "Please specify the input, skeleton and surface filenames." << endl;
		return EXIT_FAILURE;
	}

	// std::string datafile = "./region_simplification/bernardSimpleRegions.vtk";
	std::string datafile = argv[1];
	vtkNew<vtkDataSetReader> dataReader;
	dataReader->SetFileName(datafile.c_str());
	dataReader->Update();

	vtkSmartPointer<vtkDataSet> dataset = dataReader->GetOutput();

	std::string arr_name = "RegionIds";
	vtkSmartPointer<vtkDataArray> region_ids_array = dataset->GetCellData()->GetArray(arr_name.c_str());

	vtkNew<vtkIdList> unique_region_ids;
	vtkIdType numCells = region_ids_array->GetNumberOfTuples();
	for (int i = 0; i < numCells; i++)
	{
		vtkIdType id = static_cast<vtkIdType>(region_ids_array->GetTuple1(i));
		unique_region_ids->InsertUniqueId(id);
	}
	unique_region_ids->Sort();


	std::vector<Region<vtkUnstructuredGrid>> datasets;
	// First divide the ids into blocks and get the datasets into the vector
	if (unique_region_ids->GetNumberOfIds() > 5 * num_threads)
	{
		int stepSize = 5 * num_threads; // numRegions / numThreads;
		int startRegionIdx = 0, endRegionIdx = 0;
		int startRegionId = 0, endRegionId = 0;

		vtkSmartPointer<vtkThreshold> thresholdFil = vtkSmartPointer<vtkThreshold>::New();
		thresholdFil->SetInputData(dataset);
		thresholdFil->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, arr_name.c_str());

		int stepId = -1;
		bool toBreak = false;
		while (true)
		{
			stepId++;
			startRegionIdx = endRegionIdx;
			endRegionIdx = startRegionIdx + stepSize;
			startRegionId = unique_region_ids->GetId(startRegionIdx);
			endRegionId = unique_region_ids->GetId(endRegionIdx - 1);

			if (endRegionIdx >= unique_region_ids->GetNumberOfIds())
			{
				endRegionIdx = unique_region_ids->GetNumberOfIds();
				endRegionId = unique_region_ids->GetId(endRegionIdx - 1);
				toBreak = true;
			}

			cout << stepId << " " << startRegionId << " " << endRegionId << endl;

			thresholdFil->SetLowerThreshold(startRegionId);
			thresholdFil->SetUpperThreshold(endRegionId);
			thresholdFil->Update();

			// vtkSmartPointer<vtkUnstructuredGrid> tempGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
			// tempGrid->DeepCopy(thresholdFil->GetOutput());
			// cout << "Creating new region object..." << endl;

			// Store this region into a vector for later retrival
			Region<vtkUnstructuredGrid> region(stepId, thresholdFil->GetOutput());
			datasets.push_back(region);
			// cout << "Adding new region object to the list..." << endl;
			if (toBreak)
			{
				break;
			}
		}
	}
	else
	{
		vtkSmartPointer<vtkUnstructuredGrid> tempGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
		tempGrid->DeepCopy(dataset);
		Region<vtkUnstructuredGrid> region(0, tempGrid);
		datasets.push_back(region);
	}

	vtkSmartPointer<vtkPolyData> full_skeleton_data = nullptr;
	vtkSmartPointer<vtkPolyData> full_surface_data = nullptr;
	for (int datasetId = 0; datasetId < datasets.size(); datasetId++)
	{
		vtkSmartPointer<vtkPolyData> skeleton_data = nullptr;
		vtkSmartPointer<vtkPolyData> surface_data = nullptr;
		// std::string criteriaArrName = "lambda2-vtk";

		vtkUnstructuredGrid* tempGrid0 = datasets.at(datasetId).RegionData;
		vtkSmartPointer<vtkThreshold> threshold = vtkSmartPointer<vtkThreshold>::New();
		threshold->SetInputData(tempGrid0);
		threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, arr_name.c_str());

		double idsRange[2];
		tempGrid0->GetCellData()->GetArray(arr_name.c_str())->GetRange(idsRange);

		for (int id = idsRange[0]; id <= idsRange[1]; id++)
		{
			/*
			if (unique_region_ids[i] != 457)
			{
				continue;
			}
			*/
			vtkIdType idloc = unique_region_ids->FindIdLocation(id);
			if (idloc == -1)
			{
				continue;
			}
			cout << "Region Number: " << id << endl;
			threshold->SetLowerThreshold(id);
			threshold->SetUpperThreshold(id);
			threshold->Update();
			vtkSmartPointer<vtkDataSet> region_data = threshold->GetOutput();
			vtkSmartPointer<vtkPolyData> skeleton = vtkSmartPointer<vtkPolyData>::New();
			vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();
			extractSkeleton(region_data, skeleton, surface, id);

			if (skeleton == nullptr || surface == nullptr)
			{
				continue;
			}

			if (skeleton_data == nullptr)
			{
				skeleton_data = vtkSmartPointer<vtkPolyData>::New();
				skeleton_data->DeepCopy(skeleton);
				surface_data = vtkSmartPointer<vtkPolyData>::New();
				surface_data->DeepCopy(surface);
			}
			else
			{
				vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
				appendFilter->AddInputData(skeleton_data);
				appendFilter->AddInputData(skeleton);
				appendFilter->Update();
				skeleton_data->DeepCopy(appendFilter->GetOutput());

				appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
				appendFilter->AddInputData(surface_data);
				appendFilter->AddInputData(surface);
				appendFilter->Update();
				surface_data->DeepCopy(appendFilter->GetOutput());
			}
		}

		if (full_skeleton_data == nullptr)
		{
			full_skeleton_data = vtkSmartPointer<vtkPolyData>::New();
			full_skeleton_data->DeepCopy(skeleton_data);
			full_surface_data = vtkSmartPointer<vtkPolyData>::New();
			full_surface_data->DeepCopy(surface_data);
		}
		else
		{
			vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
			appendFilter->AddInputData(full_skeleton_data);
			appendFilter->AddInputData(skeleton_data);
			appendFilter->Update();
			full_skeleton_data->DeepCopy(appendFilter->GetOutput());

			appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
			appendFilter->AddInputData(full_surface_data);
			appendFilter->AddInputData(surface_data);
			appendFilter->Update();
			full_surface_data->DeepCopy(appendFilter->GetOutput());
		}

	}

	vtkNew<vtkCleanPolyData> clean;
	clean->SetInputData(full_skeleton_data);
	clean->Update();

	vtkNew<vtkDataSetWriter> writer;
	writer->SetInputData(clean->GetOutput());
	// std::string temp = "./skeleton_extraction/bernardSkeleton.vtk";
	std::string temp = argv[2];
	writer->SetFileName(temp.c_str());
	writer->SetFileTypeToBinary();
	writer->Write();

	clean->SetInputData(full_surface_data);
	clean->Update();
	writer->SetInputData(clean->GetOutput());
	// temp = "./skeleton_extraction/bernardSurface.vtk";
	temp = argv[3];
	writer->SetFileName(temp.c_str());
	writer->SetFileTypeToBinary();
	writer->Write();

	return EXIT_SUCCESS;
}