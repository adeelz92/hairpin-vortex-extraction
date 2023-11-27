// In this version of checkpoint, we use connectivity filter to grow regions using the
// local max seed points.
// Uses version 2 of get_histograms.py
// Run command
// Region_Growing.exe "Path\To\FullDataVTKfile.vtk" "Path\To\Dataset_steps.csv" "Path\To\Dataset_Regions.vtk"
#include <vtkNew.h>
#include <vtkDataSet.h>
#include <vtkOutlineFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkCubeAxesActor.h>
#include <vtkProperty.h>
#include <vtkCamera.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSMPTools.h>
#include <vtkDataSetReader.h>
#include <vtkNamedColors.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkTextProperty.h>
#include <vtkThreshold.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include "myConnectivityFilter.h"
#include <vtkAppendFilter.h>
#include <vtkDataSetWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkStaticCellLocator.h>

#include <chrono>
#include <conio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>


struct MyRange
{
	double min;
	double max;


	MyRange()
	{
		min = 0;
		max = 0;
	}

	MyRange(double _min, double _max)
	{
		min = _min;
		max = _max;
	}
};

void getCellsList(vtkDataSet* dataset, std::vector<std::array<double, 3>> positions, vtkStaticCellLocator* CellLocator,
	vtkIdList* ptsList, vtkDataArray* lambda2, double scalar_threshold)
{
	size_t size = positions.size();
	for (int pointId = 0; pointId < size; pointId++)
	{
		double *position = positions.at(pointId).data();
		// Get cell id
		vtkIdType Id = CellLocator->FindCell(position);
		if (Id >= 0)
		{
			// Get cell
			vtkNew<vtkGenericCell> cell;
			dataset->GetCell(Id, cell);

			// Get Cell Points
			vtkIdList* cellPtIds = cell->GetPointIds();

			for (int k = 0; k < cellPtIds->GetNumberOfIds(); k++)
			{
				vtkIdType id = cellPtIds->GetId(k);
				double value = *lambda2->GetTuple(id);

				if (value <= scalar_threshold)
				{
					ptsList->InsertUniqueId(Id);
					break;
				}
			}
		}
	}
}


double getAvgLambda2(vtkDataArray* lambda2)
{
	double sumLambda2 = 0;
	int totalTuples = lambda2->GetNumberOfTuples();
	for (vtkIdType tupleId = 0; tupleId < totalTuples; ++tupleId)
	{
		double value = *lambda2->GetTuple(tupleId);

		if (value < 0)
		{
			sumLambda2 += value;
		}

	}

	return (sumLambda2 / totalTuples);
}


float getMaxLambda2(vtkDataArray* lambda2)
{
	float MaxLambda2 = -INFINITY;
	for (vtkIdType tupleId = 0; tupleId < lambda2->GetNumberOfTuples(); ++tupleId)
	{
		float thisLambda2 = *lambda2->GetTuple(tupleId);
		if (thisLambda2 > MaxLambda2)
		{
			MaxLambda2 = thisLambda2;
		}
	}

	return MaxLambda2;
}


float getMinLambda2(vtkDataArray* lambda2)
{
	float MinLambda2 = INFINITY;
	for (vtkIdType tupleId = 0; tupleId < lambda2->GetNumberOfTuples(); ++tupleId)
	{
		float thisLambda2 = *lambda2->GetTuple(tupleId);
		if (thisLambda2 < MinLambda2)
		{
			MinLambda2 = thisLambda2;
		}

	}

	return MinLambda2;
}


float getMedLambda2(vtkDataArray* lambda2)
{
	float MinLambda2 = INFINITY;
	float MaxLambda2 = -INFINITY;
	for (vtkIdType tupleId = 0; tupleId < lambda2->GetNumberOfTuples(); ++tupleId)
	{
		float thisLambda2 = *lambda2->GetTuple(tupleId);
		if (thisLambda2 < MinLambda2)
		{
			MinLambda2 = thisLambda2;
		}
		if (thisLambda2 > MaxLambda2)
		{
			MaxLambda2 = thisLambda2;
		}
	}

	return (MinLambda2 + MaxLambda2) / 2;
}


std::vector<double> linspace(double start_in, double end_in, int num_in)
{
	std::vector<double> linspaced;

	if (num_in == 0)
	{
		return linspaced;
	}

	if (num_in == 1)
	{
		linspaced.push_back(start_in);
	}

	double delta = (end_in - start_in) / (num_in - 1);

	for (int i = 0; i < num_in - 1; ++i)
	{
		linspaced.push_back(start_in + delta * i);
	}
	linspaced.push_back(end_in);

	return linspaced;
}


void filterSmallRegions(vtkDataSet* input, int numRegions, int cellsthreshold,
	vtkUnstructuredGrid* output, int numThreads = 1)
{
	vtkNew<vtkAppendFilter> appendFilter;
	vtkSMPTools::LocalScope(vtkSMPTools::Config{ numThreads, "STDThread", true },
		[&]() {
		vtkSMPTools::For(0, numRegions, [&](vtkIdType regionId, vtkIdType endId) {
			vtkSMPThreadLocal<vtkSmartPointer<vtkUnstructuredGrid>> L;
			auto& dataset = L.Local();
			dataset = vtkSmartPointer<vtkUnstructuredGrid>::New();
			dataset->DeepCopy(input);
			for (; regionId < endId; ++regionId) {

				vtkNew<vtkThreshold> threshold;
				threshold->SetInputData(dataset);
				threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "RegionId");
				threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "RegionId");
				threshold->SetLowerThreshold(regionId);
				threshold->SetUpperThreshold(regionId);
				threshold->Update();

				cout << regionId << " " << threshold->GetOutput()->GetNumberOfCells() << " " << cellsthreshold << endl;

				if (threshold->GetOutput()->GetNumberOfCells() >= cellsthreshold)
				{
					vtkNew<vtkUnstructuredGrid> rData;
					rData->DeepCopy(threshold->GetOutput());

					appendFilter->AddInputData(rData);
					// appendFilter->Update();
				}
			}
		});
	});


	appendFilter->Update();
	output->DeepCopy(appendFilter->GetOutput());
}


int getDims(vtkDataSet* dataset, double bounds[6], const char coor)
{
	int count = 0;
	int numPts = dataset->GetNumberOfPoints();

	if (coor == 'x')
	{
		for (int id = 0; id < numPts; id++)
		{
			double pt[3];
			dataset->GetPoint(id, pt);

			if (pt[1] == 0 && pt[2] == 0)
			{
				count++;
			}
		}
	}

	if (coor == 'y')
	{
		for (int id = 0; id < numPts; id++)
		{
			double pt[3];
			dataset->GetPoint(id, pt);

			if (pt[0] == 0 && pt[2] == 0)
			{
				count++;
			}
		}
	}

	if (coor == 'z')
	{
		for (int id = 0; id < numPts; id++)
		{
			double pt[3];
			dataset->GetPoint(id, pt);

			if (pt[0] == 0 && pt[1] == 0)
			{
				count++;
			}
		}
	}

	return count;
}


void getCellsList(int window_size, int xDims, int yDims, int zDims,
	vtkDataSet* input, double scalar_threshold, vtkIdList * cellsList)
{
	vtkDataArray* lambda2Array = input->GetPointData()->GetArray("lambda2-vtk");

	for (int z_idx = 0; z_idx < (zDims - (window_size - 1)); z_idx += window_size)
	{
		vtkIdType* z_pts = new vtkIdType[window_size];
		for (int ii = 0; ii < window_size; ii++)
		{
			z_pts[ii] = z_idx + ii;
		}

		for (int y_idx = 0; y_idx < (yDims - (window_size - 1)); y_idx += window_size)
		{
			vtkIdType* y_pts = new vtkIdType[window_size];
			for (int ii = 0; ii < window_size; ii++)
			{
				y_pts[ii] = y_idx + ii;
			}

			for (int x_idx = 0; x_idx < (xDims - (window_size - 1)); x_idx += window_size)
			{
				vtkIdType* x_pts = new vtkIdType[window_size];
				for (int ii = 0; ii < window_size; ii++)
				{
					x_pts[ii] = x_idx + ii;
				}

				double minValue = scalar_threshold;
				vtkIdType minId = -1;
				vtkIdType minCellId = -1;

				for (int k = 0; k < window_size; k++)
				{
					for (int j = 0; j < window_size; j++)
					{
						for (int i = 0; i < window_size; i++)
						{
							vtkIdType pointId = x_pts[i] + (xDims * y_pts[j]) + (xDims * yDims * z_pts[k]);
							double value = lambda2Array->GetTuple1(pointId);

							if (value <= minValue)
							{
								vtkNew<vtkIdList> ptCellIds;
								input->GetPointCells(pointId, ptCellIds);

								minValue = value;
								minId = pointId;
								minCellId = ptCellIds->GetId(0);

								// cout << minValue << " " << minId << " " << minCellId << endl;
							}
						}
					}
				}
				delete[] x_pts;

				if (minCellId != -1)
				{
					cellsList->InsertUniqueId(minCellId);
				}
			}
			delete[] y_pts;
		}
		delete[] z_pts;
	}
}


std::vector<double> readLambda2Steps(std::string filepath, double threshold = 90)
{
	fstream fin;
	fin.open(filepath.c_str(), std::ios::in);

	if (!fin.is_open()) throw std::runtime_error("Could not open file");

	std::string line, word;
	std::vector<double> lambda2Values;
	std::vector<std::string> row;
	double val;

	while (std::getline(fin, line))
	{
		std::stringstream ss(line);

		std::vector<double> results;

		while (ss >> val) {

			results.push_back(val);

			// If the next token is a comma, ignore it and move on
			if (ss.peek() == ',') ss.ignore();
		}

		double index = results.at(0), lambda2Value = results.at(2), count = results.at(3);

		lambda2Values.push_back(lambda2Value);

	}

	return lambda2Values;
}


int main(int argc, char* argv[])
{
	vtkSMPTools::SetBackend("STDThread");
	// vtkSMPTools::Initialize(2);
	int num_threads = vtkSMPTools::GetEstimatedNumberOfThreads();
	cout << "Total threads: " << num_threads << endl;
	float cut_off;
	if (argc == 5)
	{
		cut_off = std::stof(argv[4]) / 100;
	}
	else
	{
		cut_off = 0.01 / 100;
	}

	if (argc < 4)
	{
		std::cerr << "Please specify input, csv and the output file." << endl;
		exit(1);
	}
	// std::string datafile = "./fort_data/fort180_full_float_bin.vtk";
	std::string datafile = argv[1];
	vtkNew<vtkDataSetReader> dataReader;
	dataReader->SetFileName(datafile.c_str());
	dataReader->Update();

	vtkNew<vtkTransform> transform;
	transform->Scale(1., 1., 1.);
	vtkNew<vtkTransformFilter> transformModel;
	transformModel->SetTransform(transform);
	transformModel->SetInputConnection(dataReader->GetOutputPort());
	transformModel->Update();

	vtkSmartPointer<vtkStructuredGrid> dataset = vtkStructuredGrid::SafeDownCast(transformModel->GetOutput());
	dataset->GetPointData()->SetActiveScalars("lambda2-vtk");

	double scalar_range[2];
	dataset->GetScalarRange(scalar_range);

	vtkDataArray* lambda2 = dataset->GetPointData()->GetArray("lambda2-vtk");
	// double avgLambda2 = getAvgLambda2(lambda2);
	// cout << "Avg lambda2: " << avgLambda2 << endl;

	// std::string stepsFile = "./initial_value_problem/fortFull_bin.csv";
	std::string stepsFile = argv[2];
	std::vector<double> lambda2Values = readLambda2Steps(stepsFile, 90);

	scalar_range[1] = lambda2Values.back();
	// scalar_range[1] = avgLambda2;
	cout << scalar_range[0] << " " << scalar_range[1] << endl;
	vtkNew<vtkNamedColors> colors;

	vtkNew<vtkOutlineFilter> outline;
	outline->SetInputData(dataset);
	outline->Update();
	vtkNew<vtkPolyDataMapper> outlineMapper;
	outlineMapper->SetInputConnection(outline->GetOutputPort());
	vtkNew<vtkActor> outlineActor;
	outlineActor->SetMapper(outlineMapper);
	outlineActor->GetProperty()->SetColor(colors->GetColor3d("Grey").GetData());
	outlineActor->GetProperty()->SetLineWidth(2.0);

	vtkColor3d axesActorColor = colors->GetColor3d("Grey");
	vtkNew<vtkCubeAxesActor> cubeAxesActor;
	cubeAxesActor->SetUseTextActor3D(0);
	cubeAxesActor->SetBounds(outline->GetOutput()->GetBounds());
	cubeAxesActor->SetLabelOffset(10);
	cubeAxesActor->SetTitleOffset(10);
	cubeAxesActor->SetLabelScaling(false, 0, 0, 0);
	cubeAxesActor->GetTitleTextProperty(0)->SetColor(axesActorColor.GetData());
	cubeAxesActor->GetLabelTextProperty(0)->SetColor(axesActorColor.GetData());
	cubeAxesActor->GetTitleTextProperty(1)->SetColor(axesActorColor.GetData());
	cubeAxesActor->GetLabelTextProperty(1)->SetColor(axesActorColor.GetData());
	cubeAxesActor->GetTitleTextProperty(2)->SetColor(axesActorColor.GetData());
	cubeAxesActor->GetLabelTextProperty(2)->SetColor(axesActorColor.GetData());
	cubeAxesActor->XAxisMinorTickVisibilityOff();
	cubeAxesActor->YAxisMinorTickVisibilityOff();
	cubeAxesActor->ZAxisMinorTickVisibilityOff();
	// cubeAxesActor->SetScreenSize(10.0);
	cubeAxesActor->SetFlyMode(vtkCubeAxesActor::VTK_FLY_OUTER_EDGES);

	vtkNew<vtkRenderer> ren;
	cubeAxesActor->SetCamera(ren->GetActiveCamera());
	ren->AddActor(outlineActor);
	ren->AddActor(cubeAxesActor);
	// ren->AddActor(tubeActor);
	ren->SetBackground(colors->GetColor3d("White").GetData());
	ren->GetActiveCamera()->SetPosition(0, -1, 0);
	ren->GetActiveCamera()->SetViewUp(0, 0, 1);
	ren->ResetCamera();
	ren->GetActiveCamera()->Azimuth(-30);
	ren->GetActiveCamera()->Elevation(15);
	ren->GetActiveCamera()->Zoom(1.5);

	vtkNew<vtkRenderWindow> ren_win;
	ren_win->SetSize(1000, 500);
	ren_win->SetWindowName("Vortex Core w. Q, delta and lambda Criteria");
	ren_win->AddRenderer(ren);

	vtkNew<vtkInteractorStyleTrackballCamera> irenStyle;
	vtkNew<vtkRenderWindowInteractor> iren;
	iren->SetInteractorStyle(irenStyle);
	iren->SetRenderWindow(ren_win);

	double dataBounds[6];
	dataset->GetBounds(dataBounds);
	double xMin = dataBounds[0], xMax = dataBounds[1];
	double yMin = dataBounds[2], yMax = dataBounds[3];
	double zMin = dataBounds[4], zMax = dataBounds[5];

	cout << xMin << " " << xMax << " " << yMin << " " << yMax << " " << zMin << " " << zMax << endl;
	int dims[3];
	dataset->GetDimensions(dims);
	int xDim = dims[0]; // getDims(dataset, dataBounds, 'x');
	int yDim = dims[1];  // getDims(dataset, dataBounds, 'y');
	int zDim = dims[2]; // getDims(dataset, dataBounds, 'z');

	cout << "xDim: " << xDim << endl;
	cout << "yDim: " << yDim << endl;
	cout << "zDim: " << zDim << endl;

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	vtkNew<vtkIdList> cellsList;
	getCellsList(3, xDim, yDim, zDim, dataset, scalar_range[1], cellsList);

	vtkNew<myConnectivityFilter> connectivityFilter;
	connectivityFilter->SetInputData(dataset);
	connectivityFilter->ScalarConnectivityOn();
	// connectivityFilter->FullScalarConnectivityOn();
	connectivityFilter->SetScalarRange(scalar_range);
	connectivityFilter->SetExtractionModeToCellSeededRegions();

	connectivityFilter->InitializeSeedList();
	for (vtkIdType cellIdx = 0; cellIdx < cellsList->GetNumberOfIds(); cellIdx++)
	{
		vtkIdType id = cellsList->GetId(cellIdx);
		// cout << id << " ";
		connectivityFilter->AddSeed(id);
	}
	connectivityFilter->Update();

	// int cellsthreshold = dataset->GetNumberOfCells() * 0.00001;
	int cellsthreshold = dataset->GetNumberOfCells() * cut_off;

	// Remove regions with cells less than 0.01% of total number of cells
	vtkNew<myConnectivityFilter> connectivityFilter1;
	connectivityFilter1->SetInputConnection(connectivityFilter->GetOutputPort());
	connectivityFilter1->SetExtractionModeToAllRegions();
	connectivityFilter1->ColorRegionsOn();
	connectivityFilter1->SetMinimumRegionSize(cellsthreshold);
	connectivityFilter1->Update();

	/*
	int numRegions = connectivityFilter1->GetNumberOfExtractedRegions();
	cout << "numRegions: " << numRegions << endl;

	vtkNew<vtkUnstructuredGrid> finalGrid;
	vtkSmartPointer<vtkDataSet> input = vtkDataSet::SafeDownCast(connectivityFilter1->GetOutput());
	filterSmallRegions(input, numRegions, cellsthreshold, finalGrid, 1);
	*/

	vtkNew<vtkThreshold> threshold;
	threshold->SetInputConnection(connectivityFilter1->GetOutputPort());
	threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "RegionId");
	threshold->SetLowerThreshold(1);
	threshold->SetUpperThreshold(1);
	threshold->Update();

	vtkSmartPointer<vtkDataSet> finalGrid = vtkDataSet::SafeDownCast(threshold->GetOutput());

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

	// cout << "numcells: " << finalGrid->GetNumberOfCells() << endl;

	std::string outFileName = argv[3];
	vtkNew<vtkDataSetWriter> writer;
	writer->SetInputData(finalGrid);
	writer->SetFileName(outFileName.c_str());
	writer->SetFileTypeToBinary();
	writer->Write();

	vtkNew<vtkDataSetMapper> cMapper;
	cMapper->ScalarVisibilityOff();
	cMapper->SetInputData(finalGrid);
	vtkNew<vtkActor> cActor;
	cActor->GetProperty()->SetColor(colors->GetColor3d("Orange").GetData());
	cActor->GetProperty()->SetOpacity(0.075);
	cActor->SetMapper(cMapper);
	ren->AddActor(cActor);

	ren_win->Render();
	iren->Start();

	return EXIT_SUCCESS;
}