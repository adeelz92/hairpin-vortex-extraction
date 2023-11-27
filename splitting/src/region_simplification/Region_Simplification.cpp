// In this version, we use both geometric and physical criteria to simplify regions
// In this version, we remove the unnecessary cell data after processing the regions.
// Run command
// Region_Simplification.exe "Path\To\Dataset_SplitRegions.vtk" "Path\To\Dataset_SimpleRegions.vtk"
#include <vtkActor.h>
#include <vtkAppendFilter.h>
#include <vtkCamera.h>
#include <vtkCleanPolyData.h>
#include <vtkCellData.h>
#include <vtkConnectivityFilter.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>
#include <vtkIntArray.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPointDataToCellData.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSMPTools.h>
#include <vtkStreamTracer.h>
#include <vtkThreshold.h>
#include <vtkUnstructuredGrid.h>

#include <chrono>
#include <conio.h>
#include <iostream>
#include <fstream>
#include <random>

struct Region
{
	int RegionId;
	vtkSmartPointer<vtkUnstructuredGrid> RegionData;

	Region(int regionId, vtkUnstructuredGrid* regionData)
	{
		this->RegionId = regionId;
		this->RegionData = vtkSmartPointer<vtkUnstructuredGrid>::New();
		this->RegionData->DeepCopy(regionData);
	}

	void updateData(vtkUnstructuredGrid* regionData)
	{
		this->RegionData->DeepCopy(regionData);
	}
};

std::vector<int> linspace(int start_in, int end_in, int num_in)
{
	std::vector<int> linspaced;

	if (num_in == 0)
	{
		return linspaced;
	}

	if (num_in == 1)
	{
		linspaced.push_back(start_in);
	}

	int delta = (end_in - start_in) / (num_in - 1);
	int prevId = -1;
	for (int i = 0; i < num_in - 1; ++i)
	{
		int currId = start_in + delta * i;
		if (prevId != currId)
		{
			linspaced.push_back(currId);
		}
		prevId = currId;
	}
	linspaced.push_back(end_in);

	return linspaced;
}

void getSeeds(vtkPolyData* psource, vtkDataSet* regDataset, int th = 10)
{
	/*
	// Use points as the seeds
	int totalPts = regDataset->GetNumberOfPoints();
	std::vector<int> seedIds = linspace(0, totalPts-1, int((totalPts * th) / 100));
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	double pt[3];
	for (int i = 0; i < seedIds.size(); i++)
	{
		int seedId = seedIds.at(i);
		regDataset->GetPoint(seedId, pt);
		points->InsertNextPoint(pt);
	}
	psource->SetPoints(points);
	*/

	// Use cell centers as the seeds
	int totalCells = regDataset->GetNumberOfCells();
	int numSeeds = std::max(int((totalCells * th) / 100), std::min(10, totalCells));
	std::vector<int> seedIds = linspace(0, totalCells - 1, numSeeds);
	cout << "Total Cells: " << totalCells << " Seeds Cells: " << seedIds.size() << endl;
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkNew<vtkGenericCell> cell;
	int subId;
	double pcoords[3], cellCenter[3];
	double weights[VTK_CELL_SIZE];
	for (int i = 0; i < seedIds.size(); i++)
	{
		int cellId = seedIds.at(i);
		regDataset->GetCell(cellId, cell);
		subId = cell->GetParametricCenter(pcoords);
		cell->EvaluateLocation(subId, pcoords, cellCenter, weights);
		points->InsertNextPoint(cellCenter);
	}
	psource->SetPoints(points);
}

double getCellLength(vtkCell* cell, vtkPolyData* cstmVrtxCor)
{
	vtkSmartPointer<vtkIdList> pointIds = cell->GetPointIds();
	double len = 0;
	for (vtkIdType ptIdx = 0; ptIdx < (cell->GetNumberOfPoints() - 1); ptIdx++)
	{
		vtkIdType pointId1 = pointIds->GetId(ptIdx);
		vtkIdType pointId2 = pointIds->GetId(ptIdx + 1);
		double pt1[3], pt2[3];
		cstmVrtxCor->GetPoint(pointId1, pt1);
		cstmVrtxCor->GetPoint(pointId2, pt2);
		len += vtkMath::Distance2BetweenPoints(pt1, pt2);
	}
	return len;
}

void filterSmallLines(vtkPolyData* cstmVrtxCor, vtkPolyData* output, double len_threshold)
{
	for (vtkIdType lineId = 0; lineId < cstmVrtxCor->GetNumberOfCells(); lineId++)
	{
		vtkSmartPointer<vtkCell> cell = cstmVrtxCor->GetCell(lineId);
		// double cellLen = getCellLength(cell, cstmVrtxCor);
		// cout << cell->GetLength2() << " " << len_threshold << endl;
		if (cell->GetLength2() < len_threshold)
		{
			cstmVrtxCor->DeleteCell(lineId);
		}
	}
	cstmVrtxCor->RemoveDeletedCells();

	vtkNew<vtkCleanPolyData> cleanFilter;
	cleanFilter->SetInputData(cstmVrtxCor);
	cleanFilter->Update();
	output->DeepCopy(cleanFilter->GetOutput());
}

void getCellAvgCrit(vtkDataArray* arr, vtkIdList* ptIds)
{

}

void shredRegion(vtkDataSet* tempDataset, vtkDataSet* output, std::string critArrName, double critTh,
	int mode, int* numRegs)
{
	vtkNew<vtkConnectivityFilter> maxRegExtractor;
	vtkSmartPointer<vtkThreshold> threshold = vtkSmartPointer<vtkThreshold>::New();
	vtkNew<vtkConnectivityFilter> connFilter1;
	// vtkSmartPointer<vtkDataSet> region_data;
	int numCom, numCells, totalCellNeighbors;
	vtkSmartPointer<vtkUnstructuredGrid> dataset = vtkSmartPointer<vtkUnstructuredGrid>::New();
	dataset->DeepCopy(vtkUnstructuredGrid::SafeDownCast(tempDataset));

	while (true)
	{
		vtkSmartPointer<vtkDataArray> criteriaArr = dataset->GetCellData()->GetArray(critArrName.c_str());
		numCom = criteriaArr->GetNumberOfComponents();

		// First identify the outer cells
		numCells = dataset->GetNumberOfCells();

		vtkNew<vtkIntArray> selectionId;
		selectionId->SetNumberOfComponents(1);
		selectionId->SetNumberOfTuples(numCells);
		selectionId->SetName("shreddedCell");
		dataset->GetCellData()->AddArray(selectionId);

		vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
		vtkSmartPointer<vtkIdList> cellPtIds = vtkSmartPointer<vtkIdList>::New();
		bool shredded = true;
		double values[3];
		double value;
		double mag;
		vtkSmartPointer<vtkIdList> neighborIds;
		vtkSmartPointer<vtkIdList> neighborCellIds;
		vtkSmartPointer<vtkIdList> idList;
		// Get all those points which are in the outer layer
		for (vtkIdType cellId = 0; cellId < numCells; cellId++)
		{
			dataset->GetCell(cellId, cell);
			// vtkSmartPointer<vtkGenericCell> faceCell;
			int countConnectedFaces = 0;
			for (vtkIdType faceId = 0; faceId < cell->GetNumberOfFaces(); faceId++)
			{
				int countConnectedPts = 0;
				vtkCell* faceCell = cell->GetFace(faceId);
				vtkIdList* facePtIdList = faceCell->GetPointIds();
				neighborIds = vtkSmartPointer<vtkIdList>::New();
				for (int ptId = 0; ptId < facePtIdList->GetNumberOfIds(); ptId++)
				{
					idList = vtkSmartPointer<vtkIdList>::New();
					idList->InsertNextId(facePtIdList->GetId(ptId));
					dataset->GetCellNeighbors(cellId, idList, neighborIds);
					if (neighborIds->GetNumberOfIds() > 1)
					{
						countConnectedPts++;
					}
					// cout << faceIdList->GetId(ptId) << " ";
				}

				if (countConnectedPts == 4)
				{
					countConnectedFaces++;
				}
			}

			if (mode == 1)	// Filter based on both physics and geometry
			{
				if (numCom == 1)	// Criteria Arr is a scalar
				{
					value = criteriaArr->GetTuple1(cellId);
					mag = std::sqrt(std::pow(value, 2));	// RMS of the value
				}
				else if (numCom == 3)	// Criteria Arr is a vector
				{
					criteriaArr->GetTuple(cellId, values);
					mag = vtkMath::Norm(values);
				}

				if ((countConnectedFaces <= 1 && mag < critTh / 4) || (countConnectedFaces < 1))
				{
					// 1 means cell should be shredded
					selectionId->SetTuple1(cellId, 1);
					// cout << mag << " " << critTh << " Deleted" << endl;
				}
				else
				{
					// 0 means cell should not be shredded
					selectionId->SetTuple1(cellId, 0);
					shredded = false;
				}
			}
			else if (mode == 2)	// Filter only based on geometry
			{
				if (countConnectedFaces <= 1)
				{
					// 1 means cell should be shredded
					selectionId->SetTuple1(cellId, 1);
				}
				else
				{
					// 0 means cell should not be shredded
					selectionId->SetTuple1(cellId, 0);
					shredded = false;
				}
			}

			/*
			dataset->GetCellPoints(cellId, cellPtIds);
			neighborIds = vtkSmartPointer<vtkIdList>::New();
			neighborCellIds = vtkSmartPointer<vtkIdList>::New();
			for (int i = 0; i < cellPtIds->GetNumberOfIds(); i++)
			{
				idList = vtkSmartPointer<vtkIdList>::New();
				idList->InsertNextId(cellPtIds->GetId(i));
				dataset->GetCellNeighbors(cellId, idList, neighborIds);

				for (vtkIdType j = 0; j < neighborIds->GetNumberOfIds(); j++)
				{
					neighborCellIds->InsertUniqueId(neighborIds->GetId(j));
				}
			}
			totalCellNeighbors = neighborCellIds->GetNumberOfIds();

			if (mode == 1)	// Filter based on both physics and geometry
			{
				if (numCom == 1)	// Criteria Arr is a scalar
				{
					value = criteriaArr->GetTuple1(cellId);
					mag = std::sqrt(std::pow(value, 2));	// RMS of the value
				}
				else if (numCom == 3)	// Criteria Arr is a vector
				{
					criteriaArr->GetTuple(cellId, values);
					mag = 1.0 / 2.0 * vtkMath::Norm(values);
				}

				// if ((totalCellNeighbors <= 8 && mag < critTh) || (totalCellNeighbors <= 6))
				if ((totalCellNeighbors <= 7 && mag < critTh))
				{
					// 1 means cell should be shredded
					selectionId->SetTuple1(cellId, 1);
					// cout << mag << " " << critTh << " Deleted" << endl;
				}
				else
				{
					// 0 means cell should not be shredded
					selectionId->SetTuple1(cellId, 0);
					shredded = false;
				}
			}
			else if (mode == 2)	// Filter only based on geometry
			{
				if (totalCellNeighbors <= 9)
				{
					// 1 means cell should be shredded
					selectionId->SetTuple1(cellId, 1);
				}
				else
				{
					// 0 means cell should not be shredded
					selectionId->SetTuple1(cellId, 0);
					shredded = false;
				}
			}
			*/
		}

		threshold->SetInputData(dataset);
		threshold->SetInputArrayToProcess(0, 0, 0, 1, selectionId->GetName());
		threshold->SetLowerThreshold(0);
		threshold->SetUpperThreshold(0);
		threshold->Update();

		connFilter1->SetInputData(threshold->GetOutput());
		connFilter1->SetExtractionModeToAllRegions();
		connFilter1->ColorRegionsOn();
		connFilter1->Update();
		// cout << "Regions before: " << *numRegs;
		*numRegs = connFilter1->GetNumberOfExtractedRegions();
		// cout << " Regions after: " << *numRegs << endl;

		dataset = vtkSmartPointer<vtkUnstructuredGrid>::New();
		if (*numRegs > 1)
		{
			maxRegExtractor->SetInputData(connFilter1->GetOutput());
			maxRegExtractor->SetExtractionModeToLargestRegion();
			maxRegExtractor->Update();
			dataset->DeepCopy(vtkUnstructuredGrid::SafeDownCast(maxRegExtractor->GetOutput()));
		}
		else
		{
			dataset->DeepCopy(vtkUnstructuredGrid::SafeDownCast(connFilter1->GetOutput()));
		}

		// outputGrid->DeepCopy(region_data);
		//cout << dataset->GetNumberOfCells() << endl;

		if (numCells == threshold->GetOutput()->GetNumberOfCells())
		{
			// Did not delete anything
			output->DeepCopy(dataset);
			return;
		}

		if (shredded == false)
		{
			// shredRegion(region_data, output, critArrName, critTh, 1, numRegs);
			// dataset->DeepCopy(vtkUnstructuredGrid::SafeDownCast(region_data));
		}
		else
		{
			output->DeepCopy(dataset);
			return;
		}
	}
}


int buildVortexCoreLine(vtkDataSet* regDataset, vtkDataSet* outputGrid, std::string critArrName)
{
	vtkSmartPointer<vtkPointDataToCellData> ptToCell = vtkSmartPointer<vtkPointDataToCellData>::New();
	ptToCell->PassPointDataOn();
	ptToCell->ProcessAllArraysOff();
	ptToCell->AddPointDataArray(critArrName.c_str());
	ptToCell->SetInputData(regDataset);
	ptToCell->Update();

	vtkSmartPointer<vtkUnstructuredGrid> cellData = vtkSmartPointer<vtkUnstructuredGrid>::New();
	cellData->DeepCopy(ptToCell->GetUnstructuredGridOutput());

	/**/
	// vtkSmartPointer<vtkUnstructuredGrid> outputGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkDataArray> criteriaArr = cellData->GetCellData()->GetArray(critArrName.c_str());
	int numCom = criteriaArr->GetNumberOfComponents();

	double arraySum = 0, arrayAvg = 0, mag = 0;
	double values[3];
	int nTuples = criteriaArr->GetNumberOfTuples();
	if (numCom == 1)	// Criteria Arr is a scalar
	{
		for (vtkIdType i = 0; i < nTuples; i++)
		{
			mag = criteriaArr->GetTuple1(i);
			mag = std::sqrt(std::pow(mag, 2));	// RMS of the value
			arraySum += mag;
		}
		arrayAvg = (arraySum / nTuples);
	}
	else if (numCom == 3)	// Criteria Arr is a vector
	{
		for (vtkIdType i = 0; i < nTuples; i++)
		{
			criteriaArr->GetTuple(i, values);
			mag = vtkMath::Norm(values);
			arraySum += mag;
		}
		arrayAvg = (arraySum / nTuples);
	}

	int * numRegs = new int();
	shredRegion(cellData, outputGrid, critArrName, arrayAvg, 1, numRegs);

	/*
	if (*numRegs <= 1)
	{
		return *numRegs;
	}

	vtkSmartPointer<vtkPolyData> psource = vtkSmartPointer<vtkPolyData>::New();
	getSeeds(psource, regDataset, 50);
	cout << *regDataset << endl;
	vtkSmartPointer<vtkStreamTracer> SlTracer = vtkSmartPointer<vtkStreamTracer>::New();
	SlTracer->SetInputData(regDataset);
	SlTracer->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "vorticity-oyf-vtk");
	SlTracer->SetSourceData(psource);
	SlTracer->SetMaximumPropagation(100);
	SlTracer->SetIntegratorTypeToRungeKutta45();
	SlTracer->SetIntegrationDirectionToBoth();
	SlTracer->Update();

	vtkSmartPointer<vtkPolyData> SlData = vtkSmartPointer<vtkPolyData>::New();
	double lenTh = (cellData->GetLength() * 0.5) / 100;
	filterSmallLines(SlTracer->GetOutput(), SlData, lenTh);

	vtkNew<vtkNamedColors> colors;
	vtkNew<vtkDataSetMapper> regMapper;
	regMapper->SetInputData(regDataset);
	regMapper->ScalarVisibilityOff();
	regMapper->Update();
	vtkNew<vtkActor> regActor;
	regActor->GetProperty()->SetColor(colors->GetColor3d("SkyBlue").GetData());
	regActor->SetMapper(regMapper);
	regActor->GetProperty()->SetOpacity(0.1);

	vtkNew<vtkDataSetMapper> smallRegMapper;
	smallRegMapper->SetInputData(outputGrid);
	smallRegMapper->ScalarVisibilityOff();
	smallRegMapper->Update();
	vtkNew<vtkActor> smallRegActor;
	smallRegActor->GetProperty()->SetColor(colors->GetColor3d("Green").GetData());
	smallRegActor->GetProperty()->SetLineWidth(2.0);
	smallRegActor->GetProperty()->SetOpacity(0.2);
	smallRegActor->SetMapper(smallRegMapper);

	vtkNew<vtkDataSetMapper> SlMapper;
	SlMapper->SetInputData(SlData);
	SlMapper->ScalarVisibilityOff();
	SlMapper->Update();
	vtkNew<vtkActor> SlActor;
	SlActor->GetProperty()->SetColor(colors->GetColor3d("Red").GetData());
	SlActor->GetProperty()->SetLineWidth(1);
	SlActor->GetProperty()->SetOpacity(0.3);
	SlActor->SetMapper(SlMapper);

	vtkNew<vtkRenderer> ren;
	ren->AddActor(regActor);
	ren->AddActor(smallRegActor);
	ren->AddActor(SlActor);
	ren->SetBackground(colors->GetColor3d("White").GetData());

	vtkNew<vtkRenderWindow> ren_win;
	ren_win->SetSize(1200, 600);
	ren_win->SetWindowName("Vortex Lines");
	ren_win->AddRenderer(ren);

	vtkNew<vtkRenderWindowInteractor> iren;
	vtkNew<vtkInteractorStyleTrackballCamera> irenStyle;
	iren->SetInteractorStyle(irenStyle);
	iren->SetRenderWindow(ren_win);

	ren_win->Render();
	iren->Start();
	*/

	return *numRegs;
}

int main(int argc, char* argv[])
{
	vtkSMPTools::SetBackend("STDThread");
	// vtkSMPTools::Initialize(8);
	int num_threads = vtkSMPTools::GetEstimatedNumberOfThreads();
	cout << "Total threads: " << num_threads << endl;
	if (argc < 3)
	{
		std::cerr << "Please specify input and output file names." << endl;
		return EXIT_FAILURE;
	}
	vtkNew<vtkNamedColors> colors;
	// std::string datafile = "./vortex_profile_construction/fortRegionSplitting.vtk";
	std::string datafile = argv[1];
	vtkNew<vtkDataSetReader> dataReader;
	dataReader->SetFileName(datafile.c_str());
	dataReader->Update();

	vtkSmartPointer<vtkDataSet> dataset = dataReader->GetOutput();

	std::string arr_name = "RegionIds";
	vtkSmartPointer<vtkDataArray> region_ids_array = dataset->GetCellData()->GetArray(arr_name.c_str());
	
	/*
	std::vector<int> unique_region_ids;
	vtkIdType numCells = region_ids_array->GetNumberOfTuples();
	for (vtkIdType i = 0; i < numCells; i++)
	{
		int id = region_ids_array->GetTuple1(i);
		if (std::find(unique_region_ids.begin(), unique_region_ids.end(), id) == unique_region_ids.end() && id != 0)
		{
			unique_region_ids.push_back(id);
		}
	}
	sort(unique_region_ids.begin(), unique_region_ids.end());
	*/

	vtkNew<vtkIdList> unique_region_ids;
	vtkIdType numCells = region_ids_array->GetNumberOfTuples();
	for (int i = 0; i < numCells; i++)
	{
		vtkIdType id = static_cast<vtkIdType>(region_ids_array->GetTuple1(i));
		unique_region_ids->InsertUniqueId(id);
	}
	unique_region_ids->Sort();

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	std::string criteriaArrName = "lambda2-vtk";
	int regionsAfter = 0, splitRegs = 0;

	std::vector<Region> datasets;
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

			vtkSmartPointer<vtkUnstructuredGrid> tempGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
			tempGrid->DeepCopy(thresholdFil->GetOutput());

			// Store this region into a vector for later retrival
			Region region(stepId, tempGrid);
			datasets.push_back(region);
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
		Region region(0, tempGrid);
		datasets.push_back(region);
	}

	for (int datasetId = 0; datasetId < datasets.size(); datasetId++)
	{
		vtkSmartPointer<vtkUnstructuredGrid> outDataset = nullptr;
		vtkUnstructuredGrid* tempGrid0 = datasets.at(datasetId).RegionData;
		vtkSmartPointer<vtkThreshold> threshold = vtkSmartPointer<vtkThreshold>::New();
		threshold->SetInputData(tempGrid0);
		threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, arr_name.c_str());

		double idsRange[2];
		tempGrid0->GetCellData()->GetArray(arr_name.c_str())->GetRange(idsRange);

		for (int id = idsRange[0]; id <= idsRange[1]; id++)
		{
			vtkIdType idloc = unique_region_ids->FindIdLocation(id);
			if (idloc == -1)
			{
				continue;
			}

			cout << id << endl;
			threshold->SetLowerThreshold(id);
			threshold->SetUpperThreshold(id);
			threshold->Update();
			vtkSmartPointer<vtkDataSet> region_data = threshold->GetOutput();
			vtkSmartPointer<vtkUnstructuredGrid> outputGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
			int numRegs = buildVortexCoreLine(region_data, outputGrid, criteriaArrName);

			if (numRegs > 1)
			{
				cout << "Num regs: " << numRegs << " " << id << endl;
				splitRegs++;
			}
			if (outputGrid->GetNumberOfCells() > 0)
			{
				regionsAfter++;
			}
			else
			{
				continue;
			}
			if (outDataset == nullptr)
			{
				outDataset = vtkSmartPointer<vtkUnstructuredGrid>::New();
				outDataset->DeepCopy(outputGrid);
			}
			else
			{
				vtkSmartPointer<vtkAppendFilter> appendFilter = vtkSmartPointer<vtkAppendFilter>::New();
				appendFilter->AddInputData(outDataset);
				appendFilter->AddInputData(outputGrid);
				appendFilter->Update();
				outDataset->DeepCopy(appendFilter->GetOutput());
			}
		}
		datasets.at(datasetId).updateData(outDataset);
	}

	vtkSmartPointer<vtkUnstructuredGrid> outDataset = nullptr;
	for (int datasetId = 0; datasetId < datasets.size(); datasetId++)
	{
		vtkUnstructuredGrid* outputGrid = datasets.at(datasetId).RegionData;
		if (outDataset == nullptr)
		{
			outDataset = vtkSmartPointer<vtkUnstructuredGrid>::New();
			outDataset->DeepCopy(outputGrid);
		}
		else
		{
			vtkSmartPointer<vtkAppendFilter> appendFilter = vtkSmartPointer<vtkAppendFilter>::New();
			appendFilter->AddInputData(outDataset);
			appendFilter->AddInputData(outputGrid);
			appendFilter->Update();
			outDataset->DeepCopy(appendFilter->GetOutput());
		}
	}

	/*
	for (int i = 0; i < unique_region_ids.size(); i++) {
		cout << unique_region_ids[i] << endl;
		threshold->SetLowerThreshold(unique_region_ids[i]);
		threshold->SetUpperThreshold(unique_region_ids[i]);
		threshold->Update();
		vtkSmartPointer<vtkDataSet> region_data = threshold->GetOutput();
		vtkSmartPointer<vtkUnstructuredGrid> outputGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
		int numRegs = buildVortexCoreLine(region_data, outputGrid, criteriaArrName);

		if (numRegs > 1)
		{
			cout << "Num regs: " << numRegs << " " << unique_region_ids[i] << endl;
			splitRegs++;
		}
		if (outputGrid->GetNumberOfCells() > 0)
		{
			regionsAfter++;
		}
		else
		{
			continue;
		}
		if (outDataset == nullptr)
		{
			outDataset = vtkSmartPointer<vtkUnstructuredGrid>::New();
			outDataset->DeepCopy(outputGrid);
		}
		else
		{
			vtkSmartPointer<vtkAppendFilter> appendFilter = vtkSmartPointer<vtkAppendFilter>::New();
			appendFilter->AddInputData(outDataset);
			appendFilter->AddInputData(outputGrid);
			appendFilter->Update();
			outDataset->DeepCopy(appendFilter->GetOutput());
		}
		cout << ".";
	}
	*/

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	cout << endl << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

	cout << "Total Regions: " << unique_region_ids->GetNumberOfIds() << " After: " << regionsAfter
		<< " Got Split: " << splitRegs << endl;

	int numArr = outDataset->GetCellData()->GetNumberOfArrays();
	std::vector<std::string> arrsToRemove;
	for (int i = 0; i < numArr; i++)
	{
		const char *name = outDataset->GetCellData()->GetArrayName(i);
		std::string s(name);
		std::string q = "RegionIds";
		int found = s.find(q);
		if (found != std::string::npos)
		{
			continue;
		}
		else
		{
			arrsToRemove.push_back(s);
		}
	}

	for (std::string arrName : arrsToRemove)
	{
		outDataset->GetCellData()->RemoveArray(arrName.c_str());
	}
	outDataset->Squeeze();

	// cout << *outDataset << endl;

	vtkNew<vtkDataSetWriter> writer;
	writer->SetInputData(outDataset);
	// std::string temp = "./region_simplification/fortSimpleRegions.vtk";
	std::string temp = argv[2];
	writer->SetFileName(temp.c_str());
	writer->SetFileTypeToBinary();
	writer->Write();

	return EXIT_SUCCESS;
}