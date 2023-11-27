// In this version of checkpoint, we extract the largest region in the dataset and
// try to split the region using the number of iso-surfaces components inside it.
// We use for loop function and its a breadth first approach using threading.
// If no iso-surface was computed in the previous step, that means the region
// will have no more iso-surfaces in the next step. We check this in the main
// for loop for region extraction.
// At every step, we filter small iso-surfaces. Iso-surface is considered small
// if it has less than 50 cells. This implementation uses threading to delete 
// small isosurfaces. We get the correct number of regions by using the 
// connectivity filter after filtration.
// This saves all the region id array for every level of the tree.
// Saves only those levels of RegionIds which had some splits.
// Also it uses the initial value and step size extracted from
// the histogram.
// Use get_histogram2.py strategy. 
// We filter small iso-surface based on number of cells and length of the dataset.
// In this version, we solve the wrong region id problem using the connectivity filter.
// This version use threading for the correction of wrong region ids.
// We updated the id assignment based on the center of the cell which was based on the 
// 0th point before.
// Thresholding is improved. Instead of thresholding from the whole dataset, now we
// threshold from the sub-regions.
// This version introduces the fuzzy boundary between the regions.
// This version reduce the width of the fuzzy boundary.
// Extra loop correction if performCorrections.
// In this version, I adjusted the length threshold to avoid further splitting.
// In this version, we perform correction two times.
// In this version, we fix the problem of the profile having inf values.
// Run Command
// Region_Splitting.exe "Path\To\Dataset_Regions.vtk" "Path\To\Dataset_steps.csv" "Path\To\Dataset_RegionSplitting.vtk" "Path\To\Dataset_RegionSplitting.json"
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCleanPolyData.h>
#include <vtkCharArray.h>
#include <vtkCellData.h>
#include <vtkColorSeries.h>
#include <vtkConnectivityFilter.h>
#include <vtkContourFilter.h>
#include <vtkCubeAxesActor.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>
#include <vtkGenerateGlobalIds.h>
#include <vtkFloatArray.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkjsoncpp/json/json.h>
#include <vtkLookupTable.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkOutlineFilter.h>
#include <vtkPointData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSMPTools.h>
#include <vtkStaticPointLocator.h>
#include <vtkTextProperty.h>
#include <vtkThreshold.h>
#include <vtkUnstructuredGrid.h>
#include "myPolyDataConnectivityFilter.h"

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
};


struct UniqueIdsCount
{
	int Id;
	int Count;

	UniqueIdsCount(int id)
	{
		this->Id = id;
		this->Count = 1;
	}

	void increaseCount()
	{
		this->Count++;
	}
};

class UpdateRegionIdsFunctor
{
	vtkDataSet* RegDataset;
	vtkDataSet* IsoSurfDataset;
	vtkIdTypeArray* RegColorIds;
	vtkIdType numCells;
	vtkDataArray* GlobalCellIds;
	vtkStaticPointLocator* Locator;
	vtkDataArray* SurfColorIds;
	vtkDataArray* CurrRegColorIds;

	struct LocalData
	{
		vtkSmartPointer<vtkGenericCell> Cell;
		// vtkSmartPointer<vtkPoints> CellPoints;
		// vtkSmartPointer<vtkIdList> SurfPtIdList;
		std::array<double, 3> Position;
		std::array<double, 3> Pcoords;
		std::array<double, VTK_CELL_SIZE> Weights;
		int SubId;
	};
	vtkSMPThreadLocal<LocalData> TLData;

public:

	UpdateRegionIdsFunctor(vtkDataSet* regDataset, vtkDataSet* isoSurfDataset, vtkIdTypeArray* regColorIds
		, vtkStaticPointLocator* locator, vtkDataArray* surfColorIds)
	{
		RegDataset = regDataset;
		IsoSurfDataset = isoSurfDataset;
		RegColorIds = regColorIds;
		Locator = locator;
		SurfColorIds = surfColorIds;
		numCells = RegDataset->GetNumberOfCells();
		GlobalCellIds = RegDataset->GetCellData()->GetArray("GlobalCellIds");
		CurrRegColorIds = RegDataset->GetCellData()->GetScalars();
	}

	void Initialize()
	{
		auto& tlData = this->TLData.Local();
		tlData.Cell = vtkSmartPointer<vtkGenericCell>::New();
		// tlData.CellPoints = vtkSmartPointer<vtkPoints>::New();
		// tlData.SurfPtIdList = vtkSmartPointer<vtkIdList>::New();
	}

	void operator()(vtkIdType begin, vtkIdType end)
	{
		auto& tlData = this->TLData.Local();
		auto& cell = tlData.Cell;
		// auto& cellPoints = tlData.CellPoints;
		// auto& surfPtIdList = tlData.SurfPtIdList;
		auto& position = tlData.Position;
		auto& pcoords = tlData.Pcoords;
		auto& weights = tlData.Weights;
		auto& subId = tlData.SubId;
		auto& globalCellIds = vtk::DataArrayTupleRange<1>(this->GlobalCellIds);
		auto& surfColorIds = vtk::DataArrayTupleRange<1>(this->SurfColorIds);
		auto& regColorIds = vtk::DataArrayTupleRange<1>(this->RegColorIds);
		auto& currRegColorIds = vtk::DataArrayTupleRange<1>(this->CurrRegColorIds);

		for (vtkIdType cellIdx = begin; cellIdx < end; ++cellIdx)
		{
			vtkIdType globalCellId = static_cast<vtkIdType>(globalCellIds[cellIdx][0]);
			RegDataset->GetCell(cellIdx, cell);
			// cellPoints = cell->GetPoints();
			// cellPoints->GetPoint(0, position.data());

			subId = cell->GetParametricCenter(pcoords.data());
			cell->EvaluateLocation(subId, pcoords.data(), position.data(), weights.data());

			vtkIdType surfPtId = Locator->FindClosestPoint(position.data());
			// double dat[1] = { surfColorIds[surfPtId][0] };
			// RegColorIds->SetTuple(globalCellId, dat);
			regColorIds[globalCellId][0] = surfColorIds[surfPtId][0];
			currRegColorIds[cellIdx][0] = surfColorIds[surfPtId][0];
			/*
			surfPtIdList->Reset();
			surfPtIdList->SetNumberOfIds(cellPoints->GetNumberOfPoints());

			for (vtkIdType ptId = 0; ptId < cellPoints->GetNumberOfPoints(); ptId++)
			{
				cellPoints->GetPoint(ptId, position.data());
				vtkIdType surfPtId = Locator->FindClosestPoint(position.data());
				surfPtIdList->InsertId(ptId, surfPtId);
			}

			// Check if color of all the ids are same
			vtkIdType prevSurfColor = -1;
			vtkIdType prevId;
			bool differentIds = false;
			for (vtkIdType ptId = 0; ptId < surfPtIdList->GetNumberOfIds(); ptId++)
			{
				vtkIdType surfPtId = surfPtIdList->GetId(ptId);
				vtkIdType surfColor = static_cast<vtkIdType>(surfColorIds[surfPtId][0]);

				if (prevSurfColor == -1)
				{
					prevSurfColor = surfColor;
					continue;
				}

				if (prevSurfColor == surfColor)
				{
					prevId = surfPtId;
				}
				else
				{
					differentIds = true;
					break;
				}
			}

			// cout << numCells << " " << globalCellId << endl;
			if (differentIds)
			{
				regColorIds[globalCellId][0] = 0;
				// double nullId[1] = { 0 };
				// RegColorIds->SetTuple(globalCellId, nullId);
			}
			else
			{
				regColorIds[globalCellId][0] = surfColorIds[prevId][0];
				// double dat[1] = { surfColorIds[prevId][0] };
				// RegColorIds->SetTuple(globalCellId, dat);
			}
			*/
		}
	}

	void Reduce() {}
};

/*
void UpdateCellIds(vtkDataSet* regDataset, vtkIdTypeArray* regColorIds, vtkStaticPointLocator* locator,
	vtkDataArray* surfColorIds, vtkDataArray* currRegColorIds, vtkIdList* newIds)
{
	vtkDataArray* globalCellIds = regDataset->GetCellData()->GetArray("GlobalCellIds");
	int totalCells = regDataset->GetNumberOfCells();
	vtkIdType globalCellId = -1;
	vtkSmartPointer<vtkIdList> wave = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> wave2 = vtkSmartPointer<vtkIdList>::New();
	vtkIdType* visited = new vtkIdType[totalCells];
	for (int i = 0; i < totalCells; i++)
	{
		visited[i] = -1;
	}

	vtkSmartPointer<vtkIdList> pointIds = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
	vtkNew<vtkGenericCell> cell;
	double pcoords1[3];
	double cellCenter[3];
	double weights[VTK_CELL_SIZE];
	int subId1;
	vtkIdType surfPtId;
	vtkIdType currSurfPtId, surfaceColor;
	int numPts, numCells;
	int cellId, ptId;

	struct CompletedReg
	{
		vtkIdType RegId;
		bool Completed;

		CompletedReg()
		{
			this->RegId = -1;
			this->Completed = false;
		}
	};

	int numRegs = newIds->GetNumberOfIds();
	CompletedReg* completedRegs = new CompletedReg[numRegs];
	for (int i = 0; i < numRegs; i++)
	{
		completedRegs[i].RegId = newIds->GetId(i);
		completedRegs[i].Completed = false;
	}
	bool isCompleted;

	int regId = 0;
	vtkIdType currRegId = -1;

	for (vtkIdType id = 0; id < totalCells; id++)
	{
		if (visited[id] >= 0)
		{
			continue;
		}

		wave->InsertNextId(id);
		vtkIdType numIds;
		currRegId = newIds->GetId(regId);
		regId++;

		while ((numIds = wave->GetNumberOfIds()) > 0)
		{
			for (int i = 0; i < numIds; i++)
			{
				cellId = wave->GetId(i);
				if (visited[cellId] >= 0)
				{
					continue;
				}

				regDataset->GetCell(cellId, cell);
				subId1 = cell->GetParametricCenter(pcoords1);
				cell->EvaluateLocation(subId1, pcoords1, cellCenter, weights);

				// Locate the surface id closest to the cell center
				surfPtId = locator->FindClosestPoint(cellCenter);

				surfaceColor = surfColorIds->GetTuple1(surfPtId);
				// If the surface point of the cell is not equal to the
				// currently running surface id, don't add it to the list.
				if (surfaceColor != currRegId)
				{
					continue;
				}
				else
				{
					visited[cellId] = surfPtId;
					globalCellId = static_cast<vtkIdType>(globalCellIds->GetTuple1(cellId));
					regColorIds->SetTuple1(globalCellId, surfColorIds->GetTuple1(surfPtId));
					currRegColorIds->SetTuple1(cellId, surfColorIds->GetTuple1(surfPtId));

					regDataset->GetCellPoints(cellId, pointIds);
					numPts = pointIds->GetNumberOfIds();
					for (int j = 0; j < numPts; j++)
					{
						ptId = pointIds->GetId(j);
						regDataset->GetPointCells(ptId, cellIds);
						numCells = cellIds->GetNumberOfIds();
						for (int k = 0; k < numCells; k++)
						{
							cellId = cellIds->GetId(k);
							this->Wave2->InsertNextId(cellId);
						}
					}
				}

				for (int regId = 0; regId < numRegs; regId++)
				{
					if (completedRegs[regId].RegId == surfPtId)
					{
						if (completedRegs[regId].RegId == false)
						{
							isCompleted = false;
						}
						else
						{
							isCompleted = true;
						}
						break;
					}
				}


				if (surfPtId != currSurfPtId)
				{
					continue;
				}
			}
		}
	}
}
*/

void getIsoSurfaces(vtkDataSet* dataset, float iso_value, vtkPolyData* output)
{
	// Step-1 Extract iso-surface inside the region
	vtkNew<vtkContourFilter> cFilter;
	cFilter->ComputeNormalsOff();
	cFilter->ComputeGradientsOff();
	cFilter->GenerateTrianglesOn();
	cFilter->SetInputData(dataset);
	cFilter->SetValue(0, iso_value);
	cFilter->Update();

	while (cFilter->GetOutput()->GetPointData()->GetNumberOfArrays() > 0)
	{
		cFilter->GetOutput()->GetPointData()->RemoveArray(0);
	}

	while (cFilter->GetOutput()->GetCellData()->GetNumberOfArrays() > 0)
	{
		cFilter->GetOutput()->GetCellData()->RemoveArray(0);
	}

	output->DeepCopy(cFilter->GetOutput());
}

int filterSmallIsoSurfaces(vtkPolyData* isoSurfaceData, vtkPolyData* finalIsoSurfaceData, double megaTh, 
	double regionLength, double lengthTh, int numThreads, std::ofstream& logFile, bool skipDeleting = false)
{
	// Get number of disconnected components in the iso-surface. 
	// This will be used as the indicator for how many regions, the 
	// main region should be divided.
	vtkNew<vtkPolyDataConnectivityFilter> disconnectedIsosurfaces;
	disconnectedIsosurfaces->SetInputData(isoSurfaceData);
	disconnectedIsosurfaces->SetExtractionModeToAllRegions();
	disconnectedIsosurfaces->ColorRegionsOn();
	disconnectedIsosurfaces->Update();

	vtkSmartPointer<vtkIdTypeArray> regionSizes = disconnectedIsosurfaces->GetRegionSizes();
	vtkSmartPointer<vtkPolyData> original = disconnectedIsosurfaces->GetOutput();
	vtkSmartPointer<vtkDataArray> colorIdArray = original->GetPointData()->GetScalars();
	int numRegions = disconnectedIsosurfaces->GetNumberOfExtractedRegions();
	if (skipDeleting)
	{
		finalIsoSurfaceData->DeepCopy(original);
		return numRegions;
	}
	int numCells = original->GetNumberOfCells();
	vtkSmartPointer<vtkIdList> regionsToDelete = vtkSmartPointer<vtkIdList>::New();
	int totalCellsToDelete = 0;
	if (numRegions > 5 * numThreads)
	{
		int stepSize = 5 * numThreads; // numRegions / numThreads;
		int startRegionId = 0, endRegionId = 0;
		std::vector<Region> datasets;
		vtkSmartPointer<vtkThreshold> thresholdFil = vtkSmartPointer<vtkThreshold>::New();
		thresholdFil->SetInputData(disconnectedIsosurfaces->GetOutput());
		thresholdFil->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "RegionId");
		int stepId = -1;
		bool toBreak = false;
		while (true)
		{
			stepId++;
			startRegionId = endRegionId;
			endRegionId = startRegionId + stepSize;
			if (endRegionId >= numRegions)
			{
				endRegionId = numRegions;
				toBreak = true;
			}
			thresholdFil->SetLowerThreshold(startRegionId);
			thresholdFil->SetUpperThreshold(endRegionId-1);
			thresholdFil->Update();

			cout << numRegions << " " << startRegionId << " " << endRegionId-1 << endl;
			if (logFile.is_open())
			{
				logFile << stepId << " " << startRegionId << " " << endRegionId << endl;
			}
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

		vtkSMPTools::LocalScope(vtkSMPTools::Config{ numThreads , "Sequential", true }, [&]()
		{
			vtkSMPTools::For(0, (int)datasets.size(), [&](int startId, int endId)
			{
				for (int datasetId = startId; datasetId < endId; datasetId++)
				{
					vtkUnstructuredGrid* tempGrid = datasets.at(datasetId).RegionData;
					double regionIdRange[2];
					tempGrid->GetPointData()->GetScalars()->GetRange(regionIdRange);
					vtkSmartPointer<vtkThreshold> thresholdFilter = vtkSmartPointer<vtkThreshold>::New();
					thresholdFilter->SetInputData(tempGrid);
					thresholdFilter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "RegionId");
					cout << "datasetId: " << datasetId << " ID: " << regionIdRange[0] << " Total Regions: " << regionIdRange[1] << endl;
					for (int regionId = regionIdRange[0]; regionId <= regionIdRange[1]; regionId++)
					{
						
						thresholdFilter->SetLowerThreshold(regionId);
						thresholdFilter->SetUpperThreshold(regionId);
						thresholdFilter->Update();

						double isoSurfaceLength = thresholdFilter->GetOutput()->GetLength();
						double lengthRatio = (isoSurfaceLength / regionLength) * 100;
						int size = thresholdFilter->GetOutput()->GetNumberOfCells();

						if (size < megaTh || lengthRatio < lengthTh)
						{
							totalCellsToDelete += size;
							regionsToDelete->InsertUniqueId(regionId);
							// cout << "Deleted... " << regionsToDelete->GetNumberOfIds() << endl;
						}
					}
				}
			});
		});
	}
	else
	{
		vtkNew<vtkThreshold> thresholdFilter;
		thresholdFilter->SetInputData(disconnectedIsosurfaces->GetOutput());
		thresholdFilter->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "RegionId");
		for (int regionId = 0; regionId < numRegions; regionId++)
		{
			// cout << "RegionID: " << regionId << " Total Regions: " << numRegions << endl;
			thresholdFilter->SetLowerThreshold(regionId);
			thresholdFilter->SetUpperThreshold(regionId);
			thresholdFilter->Update();

			double isoSurfaceLength = thresholdFilter->GetOutput()->GetLength();
			double lengthRatio = (isoSurfaceLength / regionLength) * 100;
			int size = thresholdFilter->GetOutput()->GetNumberOfCells();

			// vtkIdType size = static_cast<vtkIdType>(regionSizes->GetTuple1(regionId));

			// cout << "IsoSurfaceComponent: " << regionId << " LengthRatio: " << lengthRatio << " lengthTh: " << lengthTh << endl;
			// cout << "IsoSurfaceCells: " << size << " Mega: " << megaTh << endl;

			if (size < megaTh || lengthRatio < lengthTh)
			{
				totalCellsToDelete += size;
				regionsToDelete->InsertUniqueId(regionId);
				// cout << "Deleted... " << regionsToDelete->GetNumberOfIds() << endl;
			}
		}
	}

	vtkNew<vtkCleanPolyData> clean;
	vtkNew<vtkPolyDataConnectivityFilter> afterDisconnectedIsosurfaces;
	// cout << "Total Regions to Delete: " << regionsToDelete->GetNumberOfIds() << endl;
	if (regionsToDelete->GetNumberOfIds() >= 1)
	{
		vtkSmartPointer<vtkCharArray> cellsToDelete = vtkSmartPointer<vtkCharArray>::New();
		// cellsToDelete->SetNumberOfComponents(1);
		cellsToDelete->SetNumberOfTuples(numCells);
		for (vtkIdType delId = 0; delId < numCells; delId++)
		{
			char value = '0';
			cellsToDelete->SetTypedTuple(delId, &value);
		}

		int deleteCount = 0;
		vtkSMPTools::LocalScope(vtkSMPTools::Config{ numThreads, "STDThread", true },
			[&]() {
			vtkSMPTools::For(0, numCells, [&](vtkIdType startId, vtkIdType endId) {
				auto& surfColorIds = vtk::DataArrayTupleRange<1>(colorIdArray);
				auto& deleteCells = vtk::DataArrayTupleRange<1>(cellsToDelete);
				for (vtkIdType cellId = startId; cellId < endId; cellId++)
				{
					vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
					original->GetCell(cellId, cell);
					vtkIdType ptId = cell->GetPointId(0);

					vtkIdType colorId = static_cast<vtkIdType>(surfColorIds[ptId][0]);

					vtkIdType location = regionsToDelete->FindIdLocation(colorId);
					if (location != -1)
					{
						deleteCells[cellId][0] = '1';
						deleteCount++;
					}

					if (deleteCount >= totalCellsToDelete)
					{
						break;
					}

				}
			});
		});

		for (vtkIdType idx = 0; idx < numCells; idx++)
		{
			char deleted;
			cellsToDelete->GetTypedTuple(idx, &deleted);
			if (deleted == '1')
			{
				original->DeleteCell(idx);
			}
		}

		original->RemoveDeletedCells();
		disconnectedIsosurfaces->Update();

		clean->SetInputData(original);
		clean->Update();

		afterDisconnectedIsosurfaces->SetInputData(clean->GetOutput());
		afterDisconnectedIsosurfaces->SetExtractionModeToAllRegions();
		afterDisconnectedIsosurfaces->ColorRegionsOn();
		afterDisconnectedIsosurfaces->Update();
		int numRegionsBefore = numRegions;
		numRegions = afterDisconnectedIsosurfaces->GetNumberOfExtractedRegions();
		if (numRegions > 0)
		{
			cout << "Num Regs before: " << numRegionsBefore;
			cout << " Regions to Delete: " << regionsToDelete->GetNumberOfIds();
			cout << " Num Regs after: " << numRegions << endl;

			if (logFile.is_open())
			{
				logFile << "Num Regs before: " << numRegionsBefore;
				logFile << " Regions to Delete: " << regionsToDelete->GetNumberOfIds();
				logFile << " Num Regs after: " << numRegions << endl;
			}
		}

		vtkPolyData* tempData = afterDisconnectedIsosurfaces->GetOutput();
		finalIsoSurfaceData->DeepCopy(tempData);
		/*
		vtkNew<vtkPolyDataConnectivityFilter> afterDisconnectedIsosurfaces;
		afterDisconnectedIsosurfaces->SetInputData(tempData);
		afterDisconnectedIsosurfaces->SetExtractionModeToAllRegions();
		afterDisconnectedIsosurfaces->ColorRegionsOn();
		afterDisconnectedIsosurfaces->Update();

		int afterNumRegions = afterDisconnectedIsosurfaces->GetNumberOfExtractedRegions();
		if (numRegions - regionsToDelete->GetNumberOfIds() != afterNumRegions)
		{
			cout << "*****************************************" << endl;
			cout << "Problem Here: " << endl;
			cout << "Total regions before delete: " << numRegions << endl;
			cout << "Regions to delete: " << regionsToDelete->GetNumberOfIds() << endl;
			cout << "Total regions after delete: " << afterNumRegions << endl;
			cout << "##########################################" << endl;
			exit(1);
		}
		*/
		// cout << "Came here5 " << endl;
	}
	else
	{
		finalIsoSurfaceData->DeepCopy(original);
	}

	return numRegions;
}

void performCorrections(vtkDataArray* regColorIds, vtkDataSet* dataset, vtkIdList* newIds,
	double megaThreshold, int numThreads, std::ofstream& logFile)
{
	vtkDataArray* tempRegColorIds = dataset->GetCellData()->GetScalars();
	vtkDataArray* globalCellIds = dataset->GetCellData()->GetArray("GlobalCellIds");

	bool isFullDataset = true;
	if (regColorIds->GetNumberOfTuples() != tempRegColorIds->GetNumberOfTuples())
	{
		isFullDataset = false;
		// Generate new local point ids
		vtkNew<vtkIdTypeArray> localPointIds;
		localPointIds->SetNumberOfComponents(1);
		localPointIds->SetNumberOfTuples(dataset->GetNumberOfPoints());
		localPointIds->SetName("LocalPointIds");
		// Set local point ids for the new array
		for (vtkIdType pttIdx = 0; pttIdx < dataset->GetNumberOfPoints(); pttIdx++)
		{
			localPointIds->SetTuple1(pttIdx, pttIdx);
		}
		dataset->GetPointData()->AddArray(localPointIds);

		// Generate local cell ids
		vtkNew<vtkIdTypeArray> localCellIds;
		localCellIds->SetNumberOfComponents(1);
		localCellIds->SetNumberOfTuples(dataset->GetNumberOfCells());
		localCellIds->SetName("LocalCellIds");
		// Set local cell ids for the new array
		for (vtkIdType cllIdx = 0; cllIdx < dataset->GetNumberOfCells(); cllIdx++)
		{
			localCellIds->SetTuple1(cllIdx, cllIdx);
		}
		dataset->GetCellData()->AddArray(localCellIds);
	}
	/*
	cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	if (logFile.is_open())
	{
		logFile << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
	}
	*/

	std::vector<Region> datasets;
	// First divide the ids into blocks and get the datasets into the vector
	if (newIds->GetNumberOfIds() > 5 * numThreads)
	{
		int stepSize = 5 * numThreads; // numRegions / numThreads;
		int startRegionIdx = 0, endRegionIdx = 0;
		int startRegionId = 0, endRegionId = 0;

		vtkSmartPointer<vtkThreshold> thresholdFil = vtkSmartPointer<vtkThreshold>::New();
		thresholdFil->SetInputData(dataset);
		thresholdFil->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, tempRegColorIds->GetName());
		
		int stepId = -1;
		bool toBreak = false;
		while (true)
		{
			stepId++;
			startRegionIdx = endRegionIdx;
			endRegionIdx = startRegionIdx + stepSize;
			startRegionId = newIds->GetId(startRegionIdx);
			endRegionId = newIds->GetId(endRegionIdx - 1);

			if (endRegionIdx >= newIds->GetNumberOfIds())
			{
				endRegionIdx = newIds->GetNumberOfIds();
				endRegionId = newIds->GetId(endRegionIdx - 1);
				toBreak = true;
			}

			cout << stepId << " " << startRegionId << " " << endRegionId << endl;
			if (logFile.is_open())
			{
				logFile << stepId << " " << startRegionId << " " << endRegionId << endl;
			}

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

	cout << "Total sets: " << datasets.size() << endl;

	vtkSMPTools::LocalScope(vtkSMPTools::Config{ numThreads , "Sequential", true }, [&]()
	{
		vtkSMPTools::For(0, datasets.size(), [&](int startId, int endId)
		{
			for (int datasetId = startId; datasetId < endId; datasetId++)
			{
				vtkUnstructuredGrid* tempGrid0 = datasets.at(datasetId).RegionData;

				vtkSmartPointer<vtkThreshold> threshold = vtkSmartPointer<vtkThreshold>::New();
				threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, tempRegColorIds->GetName());
				threshold->SetInputData(tempGrid0);

				vtkSmartPointer<vtkConnectivityFilter> connFilter = vtkSmartPointer<vtkConnectivityFilter>::New();
				connFilter->SetExtractionModeToAllRegions();
				connFilter->ColorRegionsOn();

				double idsRange[2];
				tempGrid0->GetCellData()->GetArray(tempRegColorIds->GetName())->GetRange(idsRange);
				for (int id = idsRange[0]; id <= idsRange[1]; id++)
				{
					vtkSmartPointer<vtkIdTypeArray> maxConnRegIds = nullptr;
					vtkIdType idloc = newIds->FindIdLocation(id);
					if (idloc == -1)
					{
						continue;
					}
					vtkIdType regionId = id; // newIds->GetId(id);

					threshold->SetLowerThreshold(regionId);
					threshold->SetUpperThreshold(regionId);
					threshold->Update();

					if (threshold->GetOutput()->GetNumberOfCells() == 0)
					{
						/*
						cout << "Did not find any region: " << regionId << endl;
						if (logFile.is_open())
						{
							logFile << "Did not find any region: " << regionId << endl;
						}
						*/
						continue;
					}

					connFilter->SetInputData(threshold->GetOutput());
					connFilter->Update();

					int numRegs = connFilter->GetNumberOfExtractedRegions();
					bool isSmall = false;	// Filter very small regions as well

					/**/
					if (connFilter->GetOutput()->GetNumberOfCells() < megaThreshold)
					{
						isSmall = true;
					}
					
					if (numRegs > 1 || isSmall)
					{
						cout << "Correcting regions: " << regionId <<
							" Total Regions: " << newIds->GetNumberOfIds() <<
							" Segments: " << numRegs << endl;
						/*
						cout << "regionId: " << regionId << " numRegs: " << numRegs << endl;
						if (logFile.is_open())
						{
							logFile << "regionId: " << regionId << " numRegs: " << numRegs << endl;
						}
						*/
						if (maxConnRegIds == nullptr)
						{
							// Generate an array for maximum number of neighboring cells
							maxConnRegIds = vtkSmartPointer<vtkIdTypeArray>::New();
							maxConnRegIds->SetNumberOfComponents(1);
							maxConnRegIds->SetNumberOfTuples(dataset->GetNumberOfCells());
							maxConnRegIds->SetName("maxConnRegIds");

							// Set local point ids for the new array
							for (vtkIdType cIdx = 0; cIdx < dataset->GetNumberOfCells(); cIdx++)
							{
								maxConnRegIds->SetTuple1(cIdx, -1);
							}
						}

						std::vector<Region> regions;
						int maxRegNo = -1;
						int maxRegCells = 0;

						if (numRegs == 1 && isSmall)
						{
							vtkSmartPointer<vtkUnstructuredGrid> tempGrid1 = vtkSmartPointer<vtkUnstructuredGrid>::New();
							tempGrid1->DeepCopy(connFilter->GetOutput());

							// Store this region into a vector for later retrival
							Region region(0, tempGrid1);
							regions.push_back(region);
						}
						else
						{
							vtkNew<vtkThreshold> threshold2;
							threshold2->SetInputData(connFilter->GetOutput());
							threshold2->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "RegionId");
							for (int regNo = 0; regNo < numRegs; regNo++)
							{
								threshold2->SetLowerThreshold(regNo);
								threshold2->SetUpperThreshold(regNo);
								threshold2->Update();
								int numCells = threshold2->GetOutput()->GetNumberOfCells();
								/*
								cout << "regNo: " << regNo << " numCells: " << numCells << endl;
								if (logFile.is_open())
								{
									logFile << "regNo: " << regNo << " numCells: " << numCells << endl;
								}
								*/
								if (isSmall == false)
								{
									if (numCells > maxRegCells)
									{
										maxRegCells = numCells;
										maxRegNo = regNo;
									}
								}

								vtkSmartPointer<vtkUnstructuredGrid> tempGrid2 = vtkSmartPointer<vtkUnstructuredGrid>::New();
								tempGrid2->DeepCopy(threshold2->GetOutput());

								// Store this region into a vector for later retrival
								Region region(regNo, tempGrid2);
								regions.push_back(region);
							}

							if (maxRegCells < megaThreshold)
							{
								maxRegNo = -1;
							}
						}

						// Find the correct ids for cells with wrong region ids
						for (int regNo = 0; regNo < numRegs; regNo++)
						{
							if (regNo == maxRegNo)
							{
								continue;
							}

							vtkDataSet* RegDataset = vtkDataSet::SafeDownCast(regions.at(regNo).RegionData);
							vtkDataArray* globalPointIds;

							if (isFullDataset)
							{
								globalPointIds = RegDataset->GetPointData()->GetArray("GlobalPointIds");
							}
							else
							{
								globalPointIds = RegDataset->GetPointData()->GetArray("LocalPointIds");
							}

							int numCells = RegDataset->GetNumberOfCells();
							int numPts = RegDataset->GetNumberOfPoints();

							// Call for thread safe
							vtkSmartPointer<vtkIdList> pointCellIdsTemp = vtkSmartPointer<vtkIdList>::New();
							RegDataset->GetPointCells(0, pointCellIdsTemp);
							vtkSmartPointer<vtkIdList> gloPointCellIdsTemp = vtkSmartPointer<vtkIdList>::New();
							dataset->GetPointCells(0, gloPointCellIdsTemp);
							vtkSmartPointer<vtkIdTypeArray> LocalMaxConnRegIds = vtkSmartPointer<vtkIdTypeArray>::New();
							LocalMaxConnRegIds->DeepCopy(maxConnRegIds);

							// ##########################################################################
							vtkSMPTools::LocalScope(vtkSMPTools::Config{ numThreads, "STDThread", true }, [&]() {
								vtkSMPTools::For(0, numPts, [&](vtkIdType startId, vtkIdType endId) {
									auto& GlobalPointIds = vtk::DataArrayTupleRange<1>(globalPointIds);
									auto& GlobalCellIds = vtk::DataArrayTupleRange<1>(globalCellIds);
									auto& RegColorIds = vtk::DataArrayTupleRange<1>(regColorIds);
									auto& MaxConnRegIds = vtk::DataArrayTupleRange<1>(LocalMaxConnRegIds);

									// Get all those points which are in the outer layer
									for (vtkIdType ptId = startId; ptId < endId; ptId++)
									{
										vtkSmartPointer<vtkIdList> pointCellIds = vtkSmartPointer<vtkIdList>::New();
										RegDataset->GetPointCells(ptId, pointCellIds);

										int numPtCells = pointCellIds->GetNumberOfIds();

										// Get the boundary cells
										if (numPtCells <= 4)
										{
											vtkIdType globalPointId = static_cast<vtkIdType>(GlobalPointIds[ptId][0]);

											vtkSmartPointer<vtkIdList> gloPointCellIds = vtkSmartPointer<vtkIdList>::New();
											dataset->GetPointCells(globalPointId, gloPointCellIds);

											for (vtkIdType j = 0; j < gloPointCellIds->GetNumberOfIds(); j++)
											{
												vtkIdType gCellId = gloPointCellIds->GetId(j);
												vtkIdType globalCellId = static_cast<vtkIdType>(GlobalCellIds[gCellId][0]);
												vtkIdType colorId = static_cast<vtkIdType>(RegColorIds[globalCellId][0]);

												if (colorId == regionId)
												{
													continue;
												}

												MaxConnRegIds[gCellId][0] = colorId;
											}

										}
									}
								});
							});
							// ##########################################################################

							std::vector<UniqueIdsCount> idCounts;
							for (int d = 0; d < LocalMaxConnRegIds->GetNumberOfTuples(); d++)
							{
								vtkIdType colorId = static_cast<vtkIdType>(LocalMaxConnRegIds->GetTuple1(d));
								bool found = false;
								for (int vecId = 0; vecId < idCounts.size(); vecId++)
								{
									UniqueIdsCount & idCount = idCounts.at(vecId);
									if (idCount.Id == colorId)
									{
										// Found the id in vector
										found = true;
										idCount.increaseCount();
									}
								}

								if (found == false)
								{
									// Id not found, insert new id
									UniqueIdsCount idCount(colorId);
									idCounts.push_back(idCount);
								}
							}

							// Find the largest neighbor
							vtkIdType majorId = -1;
							vtkIdType majorNumCells = 0;
							for (int vecId = 0; vecId < idCounts.size(); vecId++)
							{
								UniqueIdsCount idCount = idCounts.at(vecId);
								if (idCount.Id == -1)
								{
									continue;
								}
								if (idCount.Count > majorNumCells)
								{
									majorNumCells = idCount.Count;
									majorId = idCount.Id;
								}
							}

							/*
							// The disconnected region is not connected to anything.
							// Delete these cells or assign them 0 id.
							if (majorId == -1)
								continue;


							cout << "regionId: " << regionId <<
								" regionCells: " << connFilter->GetOutput()->GetNumberOfCells() <<
								" regNo: " << regNo << " numCells: " << numCells <<
								" majorId: " << majorId << " majorNumCells: " << majorNumCells << endl;

							if (logFile.is_open())
							{
								logFile << "regionId: " << regionId <<
									" regionCells: " << connFilter->GetOutput()->GetNumberOfCells() <<
									" regNo: " << regNo << " numCells: " << numCells <<
									" majorId: " << majorId << " majorNumCells: " << majorNumCells << endl;
							}
							*/

							// Perform correction based on the majority cells
							vtkDataArray* GlobalCellIds = RegDataset->GetCellData()->GetArray("GlobalCellIds");
							vtkDataArray* LocalCellIds;
							vtkDataArray* localRegColorIds;
							if (isFullDataset != true)
							{
								LocalCellIds = RegDataset->GetCellData()->GetArray("LocalCellIds");
								localRegColorIds = dataset->GetCellData()->GetScalars();
							}

							for (vtkIdType cellIdx = 0; cellIdx < numCells; cellIdx++)
							{
								vtkIdType globalCellId = static_cast<vtkIdType>(GlobalCellIds->GetTuple1(cellIdx));
								if (majorId != -1)
								{
									regColorIds->SetTuple1(globalCellId, majorId);
								}
								else
								{
									regColorIds->SetTuple1(globalCellId, 0);
								}
								if (isFullDataset != true)
								{
									vtkIdType localCellId = static_cast<vtkIdType>(LocalCellIds->GetTuple1(cellIdx));
									if (majorId != -1)
									{
										localRegColorIds->SetTuple1(localCellId, majorId);
									}
									else
									{
										localRegColorIds->SetTuple1(localCellId, 0);
									}
								}
							}
						}
						/*
						cout << "##############################" << endl;
						if (logFile.is_open())
						{
							logFile << "##############################" << endl;
						}
						*/
						regions.clear();
					}
				}
			}
		});
	});
}

void introduceFuzzyBoundary(vtkDataArray* regColorIds, vtkDataSet* dataset)
{
	vtkDataArray* tempRegColorIds = dataset->GetCellData()->GetScalars();
	vtkDataArray* globalCellIds = dataset->GetCellData()->GetArray("GlobalCellIds");

	int totalCells = dataset->GetNumberOfCells();
	// vtkSmartPointer<vtkIdList> boundaryCells = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> cellPts = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> ptCells = vtkSmartPointer<vtkIdList>::New();
	vtkIdType id, jd, kd, numPts, numCells, ptId, cellId;
	vtkIdType cellColorId, currCellColorId, globalCellId;
	bool isSame;

	for (id = 0; id < totalCells; id++)
	{
		// First get the cells which have less than 8 neighbors with equal color ids
		cellColorId = tempRegColorIds->GetTuple1(id);
		if (cellColorId == 0)
		{
			continue;
		}
		dataset->GetCellPoints(id, cellPts);
		numPts = cellPts->GetNumberOfIds();
		globalCellId = globalCellIds->GetTuple1(id);

		isSame = true;
		for (jd = 0; jd < numPts; jd++)
		{
			ptId = cellPts->GetId(jd);
			dataset->GetPointCells(ptId, ptCells);
			numCells = ptCells->GetNumberOfIds();

			for (kd = 0; kd < numCells; kd++)
			{
				cellId = ptCells->GetId(kd);
				currCellColorId = tempRegColorIds->GetTuple1(cellId);
				if (currCellColorId != cellColorId && currCellColorId != 0)
				{
					tempRegColorIds->SetTuple1(id, 0);
					regColorIds->SetTuple1(globalCellId, 0);
					isSame = false;
					break;
				}
			}

			if (isSame == false)
			{
				// This is a boundary cell
				// boundaryCells->InsertNextId(id);
				break;
			}
		}
	}

	/*
	vtkIdType globalCellId, localColor, globalColor;
	totalCells = boundaryCells->GetNumberOfIds();
	for (id = 0; id < totalCells; id++)
	{
		cellId = boundaryCells->GetId(id);
		globalCellId = globalCellIds->GetTuple1(cellId);
		tempRegColorIds->SetTuple1(cellId, 0);
		regColorIds->SetTuple1(globalCellId, 0);
	}
	*/
}

int assignNextColor(vtkDataSet* regDataset, float currIsoValue, int prevNumRegions,
	vtkIdTypeArray* regColorIds, vtkIdList* newIds, double megaThreshold,
	double datasetLength, float vsf, int numThreads, std::ofstream& logFile)
{
	double regionLength = regDataset->GetLength();
	double lengthThreshold = (regionLength / datasetLength);
	lengthThreshold = vsf / lengthThreshold;

	// Step-1 Extract iso-surface inside the region	
	vtkNew<vtkPolyData> cFilterOutput;
	getIsoSurfaces(regDataset, currIsoValue, cFilterOutput);

	// cout << " megaThreshold: " << megaThreshold << endl;
	// cout << "totalLengthOfDataset: " << regDataset->GetLength() << endl;

	if (cFilterOutput->GetNumberOfCells() <= megaThreshold)
	{
		return -1;
	}

	vtkNew<vtkPolyData> isoSurfDataset;
	int numRegions = filterSmallIsoSurfaces(cFilterOutput, isoSurfDataset, 
		megaThreshold, regionLength, lengthThreshold, numThreads, logFile);
	// cout << "Total regions: " << numRegions << endl;	
	// cout << *disconnectedIsosurfaces->GetOutput() << endl;

	if (numRegions <= 1)
	{
		return 0;
	}

	if (isoSurfDataset->GetNumberOfCells() <= megaThreshold)
	{
		return -1;
	}
	vtkDataArray* surfColorIds = isoSurfDataset->GetPointData()->GetScalars();
	// Update the color ids of the iso-surfaces depending on how many values we had previously.	
	for (int ptId = 0; ptId < isoSurfDataset->GetNumberOfPoints(); ptId++)
	{
		vtkIdType surfColor = (vtkIdType)surfColorIds->GetTuple(ptId)[0];
		// cout << "Surf color: " << surfColor << " New Color: " << surfColor + prevNumRegions << endl;	
		double newSurfColor[1] = { surfColor + prevNumRegions };
		surfColorIds->SetTuple(ptId, newSurfColor);
	}

	for (int i = 0; i < surfColorIds->GetNumberOfTuples(); i++)
	{
		vtkIdType id = static_cast<vtkIdType>(surfColorIds->GetTuple1(i));
		newIds->InsertUniqueId(id);
	}
	newIds->Sort();

	vtkNew<vtkStaticPointLocator> locator;
	locator->SetDataSet(isoSurfDataset);
	locator->BuildLocator();
	UpdateRegionIdsFunctor functor(regDataset, isoSurfDataset, regColorIds, locator, surfColorIds);
	vtkSMPTools::LocalScope(vtkSMPTools::Config{ numThreads, "STDThread", true }, [&]() {
		vtkSMPTools::For(0, regDataset->GetNumberOfCells(), functor);
	});
	regColorIds->Modified();

	performCorrections(regColorIds, regDataset, newIds, megaThreshold, numThreads, logFile);
	// introduceFuzzyBoundary(regColorIds, regDataset);
	// performCorrections(regColorIds, regDataset, newIds, megaThreshold, numThreads, logFile);
	/*
	vtkDataArray* localRegionColorIds = regDataset->GetCellData()->GetScalars();
	vtkSmartPointer<vtkIdList> uniqueIdsAfterProcessing = vtkSmartPointer<vtkIdList>::New();
	for (vtkIdType tupleId = 0; tupleId < localRegionColorIds->GetNumberOfTuples(); tupleId++)
	{
		vtkIdType colorId = localRegionColorIds->GetTuple1(tupleId);
		uniqueIdsAfterProcessing->InsertUniqueId(colorId);
	}
	for (vtkIdType uId = 0; uId < uniqueIdsAfterProcessing->GetNumberOfIds(); uId++)
	{
		cout << "regionId After: " << uniqueIdsAfterProcessing->GetId(uId) << endl;
	}
	cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
	*/
	return numRegions;
}

int buildVortexProfileVectors(std::vector<float>& vortexProfile, vtkDataSet* newRegDataset)
{
	/*
	0. Avg(lambda_2)		1. Avg(oyf)				2. Size
	3. max(x) - min(x)		4. max(y) - min(y)		5. max(z) - min(z)
	*/

	float arraySum = 0;
	float arrayAvg;
	int nCount = 1;
	int nTuples = newRegDataset->GetNumberOfPoints();
	if (nTuples == 0)
	{
		return 0;
	}
	// 0. Avg(lambda_2)
	vtkDataArray* arr = newRegDataset->GetPointData()->GetArray("lambda2-vtk");
	for (vtkIdType i = 0; i < nTuples; i++)
	{
		float value = arr->GetTuple1(i);
		if (value < 0)
		{
			// We only consider negative value for lambda2
			arraySum += value;
			nCount++;
		}
	}
	arrayAvg = (arraySum / nCount);
	vortexProfile.push_back(arrayAvg);
	/*
	// 1. Avg(oyf)
	if (newRegDataset->GetPointData()->HasArray("oyf"))
	{
		arr = newRegDataset->GetPointData()->GetArray("oyf");
		arraySum = 0;
		nCount = 1;

		for (vtkIdType i = 0; i < nTuples; i++)
		{
			float value = arr->GetTuple1(i);

			if (value <= 0)
			{
				continue;
			}

			nCount++;
			arraySum += value;
		}

		arrayAvg = (arraySum / nCount);
		vortexProfile.push_back(arrayAvg);

	}
	else
	{
		arr = newRegDataset->GetPointData()->GetArray("vorticity-vtk");
		arraySum = 0;
		nCount = 1;
		for (vtkIdType i = 0; i < nTuples; i++)
		{
			double values[3];
			arr->GetTuple(i, values);

			if (values[1] <= 0)
			{
				continue;
			}

			arraySum += values[1];
			nCount++;
		}
		arrayAvg = (arraySum / nCount);
		vortexProfile.push_back(arrayAvg);
	}

	// 2. Size
	vortexProfile.push_back(static_cast<float>(newRegDataset->GetNumberOfCells()));

	// 3. max(x) - min(x)
	// 4. max(y) - min(y)
	// 5. max(z) - min(z)
	double point[3];
	float xMin = +INFINITY, xMax = -INFINITY;
	float yMin = +INFINITY, yMax = -INFINITY;
	float zMin = +INFINITY, zMax = -INFINITY;
	float x, y, z;
	for (vtkIdType i = 0; i < nTuples; i++)
	{
		newRegDataset->GetPoint(i, point);
		x = point[0], y = point[1], z = point[2];
		if (x < xMin) { xMin = x; }
		if (x > xMax) { xMax = x; }
		if (y < yMin) { yMin = y; }
		if (y > yMax) { yMax = y; }
		if (z < zMin) { zMin = z; }
		if (z > zMax) { zMax = z; }
	}

	vortexProfile.push_back(xMax - xMin);
	vortexProfile.push_back(yMax - yMin);
	vortexProfile.push_back(zMax - zMin);
	*/
	return nTuples;
}

Json::Value buildInitialTree(vtkDataArray* regColorIds, std::vector<Region> & regions,
	vtkDataSet* dataset, double lambda2, int numThreads)
{
	int level = 1;
	int parent = 0;

	vtkNew<vtkIdList> uniqueIds;
	for (int i = 0; i < regColorIds->GetNumberOfTuples(); i++)
	{
		vtkIdType id = static_cast<vtkIdType>(regColorIds->GetTuple1(i));
		uniqueIds->InsertUniqueId(id);
	}
	uniqueIds->Sort();
	std::vector<Region> datasets;
	// First divide the ids into blocks and get the datasets into the vector
	if (uniqueIds->GetNumberOfIds() > 5 * numThreads)
	{
		int stepSize = 5 * numThreads; // numRegions / numThreads;
		int startRegionIdx = 0, endRegionIdx = 0;
		int startRegionId = 0, endRegionId = 0;

		vtkSmartPointer<vtkThreshold> thresholdFil = vtkSmartPointer<vtkThreshold>::New();
		thresholdFil->SetInputData(dataset);
		thresholdFil->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, regColorIds->GetName());

		int stepId = -1;
		bool toBreak = false;
		while (true)
		{
			stepId++;
			startRegionIdx = endRegionIdx;
			endRegionIdx = startRegionIdx + stepSize;
			startRegionId = uniqueIds->GetId(startRegionIdx);
			endRegionId = uniqueIds->GetId(endRegionIdx - 1);

			if (endRegionIdx >= uniqueIds->GetNumberOfIds())
			{
				endRegionIdx = uniqueIds->GetNumberOfIds();
				endRegionId = uniqueIds->GetId(endRegionIdx - 1);
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

	Json::Value tree;
	tree["name"] = "root";
	Json::Value children;

	vtkSMPTools::LocalScope(vtkSMPTools::Config{ numThreads , "STDThread", true }, [&]()
	{
		vtkSMPTools::For(0, (int)datasets.size(), [&](int startId, int endId)
		{
			for (int datasetId = startId; datasetId < endId; datasetId++)
			{
				vtkUnstructuredGrid* tempGrid0 = datasets.at(datasetId).RegionData;

				double idsRange[2];
				tempGrid0->GetCellData()->GetArray(regColorIds->GetName())->GetRange(idsRange);

				vtkNew<vtkThreshold> threshold;
				threshold->SetInputData(tempGrid0);
				threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, regColorIds->GetName());

				for (int id = idsRange[0]; id <= idsRange[1]; id++)
				{
					vtkSmartPointer<vtkIdTypeArray> maxConnRegIds = nullptr;
					vtkIdType idloc = uniqueIds->FindIdLocation(id);
					if (idloc == -1)
					{
						continue;
					}

					vtkIdType regionId = id; // newIds->GetId(id);
					cout << "Processing regions: " << regionId <<
						" Total Regions: " << uniqueIds->GetId(uniqueIds->GetNumberOfIds() - 1) << endl;

					// Don't add regionId = 0, it indicates fuzzy boundary
					if (regionId == 0)
					{
						continue;
					}

					threshold->SetLowerThreshold(regionId);
					threshold->SetUpperThreshold(regionId);
					threshold->Update();

					vtkSmartPointer<vtkUnstructuredGrid> tempGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
					tempGrid->DeepCopy(threshold->GetOutput());

					// Store this region into a vector for later retrival
					Region region(regionId, tempGrid);
					regions.push_back(region);

					vtkSmartPointer<vtkDataSet> newRegDataset = vtkDataSet::SafeDownCast(threshold->GetOutputDataObject(0));
					Json::Value node;
					std::string name = std::to_string(regionId) + "_" + std::to_string(level) + "_" + std::to_string(parent) +
						"_" + std::to_string(lambda2);

					std::vector<float> vortexProfile;
					int numPts = buildVortexProfileVectors(vortexProfile, newRegDataset);
					if (numPts == 0)
					{
						return NULL;
					}
					int nCriteria = vortexProfile.size();

					// Add the vortex profile to the tree node name for later retrieval
					for (int ii = 0; ii < nCriteria; ii++)
					{
						name += "_" + std::to_string(vortexProfile[ii]);
					}

					node["name"] = name;
					node["children"] = {};
					children.append(node);
				}
			}
		});
	});
	tree["children"] = children;
	return tree;
}

void updateJson(Json::Value& tree, vtkIdList* newIds, vtkDataSet* dataset, vtkDataArray* regColorIds,
	std::vector<Region> & regions, int& level, int& regionId, double lambda2, int parent = 0)
{
	int numIds = newIds->GetNumberOfIds();
	std::string regionName = std::to_string(regionId);

	Json::Value &nodes = tree["children"];
	int newLevel = level + 1;

	for (unsigned i = 0; i < nodes.size(); i++)
	{
		Json::Value &node = nodes[i];
		std::string nodeName = node["name"].asString();
		std::string token = nodeName.substr(0, nodeName.find("_"));

		if (regionName == token)
		{
			vtkNew<vtkThreshold> threshold;
			threshold->SetInputData(dataset);
			threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, regColorIds->GetName());

			// If there is only id, we just need to update the node
			if (numIds == 1)
			{
				vtkIdType newRegionId = newIds->GetId(0);
				Json::Value newNode;
				std::string newName = std::to_string(newRegionId) + "_" + std::to_string(newLevel) + "_" + std::to_string(regionId)
					+ "_" + std::to_string(lambda2);

				threshold->SetLowerThreshold(newRegionId);
				threshold->SetUpperThreshold(newRegionId);
				threshold->Update();
				/*
				vtkNew<vtkConnectivityFilter> connFilter;
				connFilter->SetInputData(threshold->GetOutput());
				connFilter->SetExtractionModeToAllRegions();
				connFilter->ColorRegionsOn();
				connFilter->Update();

				int numRegs = connFilter->GetNumberOfExtractedRegions();
				cout << "Reg Id: " << newRegionId << " Connected Regs: " << numRegs << endl;
				if (numRegs > 1)
				{
					for (int regNo = 0; regNo < numRegs; regNo++)
					{
						vtkNew<vtkThreshold> threshold2;
						threshold2->SetInputData(connFilter->GetOutput());
						threshold2->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "RegionId");
						threshold2->SetLowerThreshold(regNo);
						threshold2->SetUpperThreshold(regNo);
						threshold2->Update();
						int numCells = threshold2->GetOutput()->GetNumberOfCells();
						cout << "Reg No: " << regNo << " Num Cells: " << numCells << endl;
					}
				}
				*/

				// Store this region into a vector for later retrival
				vtkSmartPointer<vtkUnstructuredGrid> tempGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
				tempGrid->DeepCopy(threshold->GetOutput());

				Region region(newRegionId, tempGrid);
				regions.push_back(region);

				vtkSmartPointer<vtkDataSet> newRegDataset = vtkDataSet::SafeDownCast(threshold->GetOutputDataObject(0));

				// Build the vortex profile
				// const int nCriteria = 36;
				// double vortexProfile[nCriteria];
				// buildVortexProfile(vortexProfile, newRegDataset);

				std::vector<float> vortexProfile;
				int numPts = buildVortexProfileVectors(vortexProfile, newRegDataset);
				if (numPts == 0)
				{
					continue;
				}
				int nCriteria = vortexProfile.size();

				// Add the vortex profile to the tree node name for later retrieval
				for (int ii = 0; ii < nCriteria; ii++)
				{
					newName += "_" + std::to_string(vortexProfile[ii]);
				}

				node["name"] = newName;
				return;
			}
			else
			{
				Json::Value children;

				for (int j = 0; j < numIds; j++)
				{
					vtkIdType newRegionId = newIds->GetId(j);
					Json::Value newNode;
					std::string newName = std::to_string(newRegionId) + "_" + std::to_string(newLevel) +
						"_" + std::to_string(regionId) + "_" + std::to_string(lambda2);

					threshold->SetLowerThreshold(newRegionId);
					threshold->SetUpperThreshold(newRegionId);
					threshold->Update();

					/*
					vtkNew<vtkConnectivityFilter> connFilter;
					connFilter->SetInputData(threshold->GetOutput());
					connFilter->SetExtractionModeToAllRegions();
					connFilter->ColorRegionsOn();
					connFilter->Update();

					int numRegs = connFilter->GetNumberOfExtractedRegions();
					if (numRegs > 1)
					{
						for (int regNo = 0; regNo < numRegs; regNo++)
						{
							vtkNew<vtkThreshold> threshold2;
							threshold2->SetInputData(connFilter->GetOutput());
							threshold2->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "RegionId");
							threshold2->SetLowerThreshold(regNo);
							threshold2->SetUpperThreshold(regNo);
							threshold2->Update();
							int numCells = threshold2->GetOutput()->GetNumberOfCells();
							cout << "Reg No: " << regNo << " Num Cells: " << numCells << endl;
						}
					}
					*/

					// Store this region into a vector for later retrival
					vtkSmartPointer<vtkUnstructuredGrid> tempGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
					tempGrid->DeepCopy(threshold->GetOutput());

					Region region(newRegionId, tempGrid);
					regions.push_back(region);

					vtkSmartPointer<vtkDataSet> newRegDataset = vtkDataSet::SafeDownCast(threshold->GetOutputDataObject(0));

					// Build the vortex profile
					// const int nCriteria = 36;
					// double vortexProfile[nCriteria];
					// buildVortexProfile(vortexProfile, newRegDataset);

					std::vector<float> vortexProfile;
					int numPts = buildVortexProfileVectors(vortexProfile, newRegDataset);
					if (numPts == 0)
					{
						continue;
					}
					int nCriteria = vortexProfile.size();

					// Add the vortex profile to the tree node name for later retrieval
					for (int ii = 0; ii < nCriteria; ii++)
					{
						newName += "_" + std::to_string(vortexProfile[ii]);
					}

					newNode["name"] = newName;
					newNode["children"] = {};
					children.append(newNode);

				}
				node["children"] = children;
				return;
			}
		}
		else
		{
			updateJson(nodes[i], newIds, dataset, regColorIds, regions, level, regionId,
				lambda2, parent);
		}
	}
}

int colorRegions(vtkDataSet* dataset, std::vector<double> lambda2Values, float vsf, 
	int numThreads, std::ofstream& logFile, Json::Value& tree)
{
	/*
	Json::Reader json_reader;
	std::ifstream f;
	f.open("flare.json");

	Json::Value tree;
	if (!json_reader.parse(f, tree))
	{
		cout << json_reader.getFormattedErrorMessages() << endl;
	}


	updateJson(tree, "Scale");
	*/

	vtkNew<vtkIdTypeArray> regColorIds;
	regColorIds->SetNumberOfComponents(1);
	regColorIds->SetNumberOfTuples(dataset->GetNumberOfCells());
	regColorIds->SetName("RegionIds");

	// Set all region ids to 0
	for (vtkIdType cellIdx = 0; cellIdx < dataset->GetNumberOfCells(); cellIdx++)
	{
		double nullId[1] = { 0 };
		regColorIds->SetTuple(cellIdx, nullId);
		// regColorIds->SetTypedTuple(cellIdx, { 0 });
	}

	dataset->GetCellData()->AddArray(regColorIds);
	dataset->GetCellData()->SetActiveScalars(regColorIds->GetName());
	
	double totalCellsInDataset = dataset->GetNumberOfCells();
	double megaThreshold = 50;
	double datasetLength = dataset->GetLength();
	double lengthThreshold = (datasetLength / datasetLength);
	lengthThreshold = 1 / lengthThreshold;
	double first_iso_value = lambda2Values.back();
	lambda2Values.pop_back();

	vtkNew<vtkPolyData> cFilterOutput;
	getIsoSurfaces(dataset, first_iso_value, cFilterOutput);

	if (cFilterOutput->GetNumberOfPoints() < 1)
	{
		return 0;
	}

	vtkNew<vtkPolyData> isoSurfDataset;
	int numRegions = filterSmallIsoSurfaces(cFilterOutput, isoSurfDataset, megaThreshold, datasetLength,
		lengthThreshold, numThreads, logFile, true) + 1;	// Including 0

	vtkDataArray* surfColorIds = isoSurfDataset->GetPointData()->GetScalars();
	// Update the color ids of the iso-surfaces depending on how many values we had previously.
	for (int ptId = 0; ptId < isoSurfDataset->GetNumberOfPoints(); ptId++)
	{
		vtkIdType surfColor = (vtkIdType)surfColorIds->GetTuple(ptId)[0];
		double newSurfColor[1] = { surfColor + 1 };
		surfColorIds->SetTuple(ptId, newSurfColor);
	}
	
	vtkSmartPointer<vtkIdList> initialIds = vtkSmartPointer<vtkIdList>::New();
	for (int i = 0; i < surfColorIds->GetNumberOfTuples(); i++)
	{
		vtkIdType id = static_cast<vtkIdType>(surfColorIds->GetTuple1(i));
		initialIds->InsertUniqueId(id);
	}
	initialIds->Sort();

	// Build a locator to locate the closest points
	vtkNew<vtkStaticPointLocator> locator;
	locator->SetDataSet(isoSurfDataset);
	locator->BuildLocator();
	cout << "Assigning color ids..." << endl;
	UpdateRegionIdsFunctor functor(dataset, isoSurfDataset, regColorIds, locator, surfColorIds);
	vtkSMPTools::LocalScope(vtkSMPTools::Config{ numThreads, "STDThread", true }, [&]() {
		vtkSMPTools::For(0, dataset->GetNumberOfCells(), functor);
	});
	
	cout << "Finished assigning color ids..." << endl;
	performCorrections(regColorIds, dataset, initialIds, megaThreshold, numThreads, logFile);
	cout << "Finished performing corrections..." << endl;
	regColorIds->Modified();
	// introduceFuzzyBoundary(regColorIds, dataset);
	// performCorrections(regColorIds, dataset, initialIds, megaThreshold, numThreads, logFile);
	/*
	vtkDataArray* localRegionColorIds = dataset->GetCellData()->GetScalars();
	vtkSmartPointer<vtkIdList> uniqueIdsAfterProcessing = vtkSmartPointer<vtkIdList>::New();
	for (vtkIdType tupleId = 0; tupleId < localRegionColorIds->GetNumberOfTuples(); tupleId++)
	{
		vtkIdType colorId = localRegionColorIds->GetTuple1(tupleId);
		uniqueIdsAfterProcessing->InsertUniqueId(colorId);
	}
	for (vtkIdType uId = 0; uId < uniqueIdsAfterProcessing->GetNumberOfIds(); uId++)
	{
		cout << "regionId After: " << uniqueIdsAfterProcessing->GetId(uId) << endl;
	}
	cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
	
	double first_iso_value = lambda2Values.back();
	*/

	std::vector<Region> regions;
	tree = buildInitialTree(regColorIds, regions, dataset, first_iso_value, numThreads);
	cout << "Finished building initial tree.." << endl;
	if (logFile.is_open())
	{
		logFile << "Finished building initial tree.." << endl;
	}
	
	std::vector<int> checkRegions;
	double range[2];
	regColorIds->GetRange(range);
	cout << "Num Regs: " << numRegions << " Iso-Value: " << first_iso_value << " Range: " << range[0] << " " << range[1] << endl;
	if (logFile.is_open())
	{
		logFile << "Num Regs: " << numRegions << " Iso-Value: " << first_iso_value << " Range: " << range[0] << " " << range[1] << endl;
	}
	// double lambda2Range[2];
	// dataset->GetPointData()->GetScalars()->GetRange(lambda2Range);
	int level = 1;
	std::string arrayName = "RegionIds_" + std::to_string(level);

	vtkNew<vtkIdTypeArray> prevRegColorIds;
	prevRegColorIds->DeepCopy(regColorIds);
	prevRegColorIds->SetName(arrayName.c_str());
	dataset->GetCellData()->AddArray(prevRegColorIds);
	dataset->GetCellData()->SetActiveScalars(prevRegColorIds->GetName());

	while (lambda2Values.size() > 0)
	{
		double isoValue = lambda2Values.back();
		lambda2Values.pop_back();

		// cout << arrayName << endl;
		vtkNew<vtkIdList> uniqueIds;
		for (int i = 0; i < regColorIds->GetNumberOfTuples(); i++)
		{
			vtkIdType id = static_cast<vtkIdType>(regColorIds->GetTuple1(i));
			if (id == 0)
			{
				continue;
			}
			uniqueIds->InsertUniqueId(id);
		}
		uniqueIds->Sort();

		vtkDataSet* newRegDataset;
		vtkNew<vtkThreshold> threshold;
		threshold->SetInputData(dataset);
		threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, arrayName.c_str());
		// cout << *dataset->GetCellData() << endl;
		int currenRegions = numRegions;
		int totalRegionsthisIteration = 0;
		for (int i = 0; i < uniqueIds->GetNumberOfIds(); i++)
		{
			int regionId = uniqueIds->GetId(i);
			if (checkRegions.size() >= regionId)
			{
				int numPrevRegions = checkRegions.at(regionId - 1);
				if (numPrevRegions == -1)
				{
					// Reducing lambda2 value won't introduce new regions because the number 
					// of iso-surfaces are already 0

					continue;
				}
			}

			bool foundInVec = false;
			int vecId = 0;
			// Try to find the regionId in the saved vector
			for (; vecId < regions.size(); vecId++)
			{
				if (regionId == regions.at(vecId).RegionId)
				{
					foundInVec = true;
					break;
				}
			}

			/**/
			if (foundInVec)
			{
				// cout << *regions.at(vecId).RegionData << endl;
				newRegDataset = vtkDataSet::SafeDownCast(regions.at(vecId).RegionData);
			}
			else
			{
				threshold->SetLowerThreshold(regionId);
				threshold->SetUpperThreshold(regionId);
				threshold->Update();
				newRegDataset = vtkDataSet::SafeDownCast(threshold->GetOutputDataObject(0));

				// Store this region into a vector for later retrival
				// vtkSmartPointer<vtkUnstructuredGrid> tempGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
				// tempGrid->DeepCopy(threshold->GetOutput());
				Region temp(regionId, threshold->GetOutput());
				regions.push_back(temp);
				// cout << "Not found in vec: " << regionId << endl;
			}
			if (newRegDataset->GetNumberOfCells() == 0)
			{
				continue;
			}
			vtkSmartPointer<vtkIdList> newIds = vtkSmartPointer<vtkIdList>::New();
			int newRegions = assignNextColor(newRegDataset, isoValue, numRegions, regColorIds, newIds, megaThreshold,
				datasetLength, vsf, numThreads, logFile);
			if (checkRegions.size() < regionId)
			{
				checkRegions.push_back(newRegions);
			}
			else
			{
				checkRegions.at(regionId - 1) = newRegions;
			}

			if (newRegions <= 0)
			{
				newRegions = 0;
				continue;
			}

			if (newIds->GetNumberOfIds() > 0)
			{
				vtkDataArray* localRegColorIds = newRegDataset->GetCellData()->GetScalars();
				updateJson(tree, newIds, newRegDataset, localRegColorIds, regions, level, regionId, isoValue);
				/**/
				if (foundInVec)
				{
					// cout << "Found in Vec: " << vecId << " " << regionId << " " << newRegions << " " << newIds->GetNumberOfIds() << endl;
					regions.erase(regions.begin() + vecId);
					// regions.at(vecId).RegionData->Delete();
					// regions.at(vecId).RegionId = -1;
				}
				
			}

			totalRegionsthisIteration += newRegions;
			regColorIds->GetRange(range);

			/**/
			cout << regionId << " Old Regs: " << numRegions << " New Regs: " << newRegions << " Total Regs: "
				<< numRegions + newRegions << " Iso: " << isoValue << " Range: " << range[0] << " " << range[1] << endl;
			
			if (logFile.is_open())
			{
				logFile << regionId << " Old Regs: " << numRegions << " New Regs: " << newRegions << " Total Regs: "
					<< numRegions + newRegions << " Iso: " << isoValue << " Range: " << range[0] << " " << range[1] << endl;
			}
			numRegions += newRegions;
			// break;
		}

		cout << "***************************************************************" << endl;
		if (logFile.is_open())
		{
			logFile << "***************************************************************" << endl;
		}
		// Update the array name if some region actually got split otherwise there is
		// no need to add a new array
		if (totalRegionsthisIteration > 0)
		{
			level += 1;
			arrayName = "RegionIds_" + std::to_string(level);
			vtkNew<vtkIdTypeArray> newRegColorIds;
			newRegColorIds->DeepCopy(regColorIds);
			newRegColorIds->SetName(arrayName.c_str());
			dataset->GetCellData()->AddArray(newRegColorIds);
			dataset->GetCellData()->SetActiveScalars(newRegColorIds->GetName());
		}

		// break;
		/*
		vtkNew<vtkDataSetWriter> writer;
		writer->SetInputData(dataset);
		std::string temp = "./selectData/regionSplit";
		std::string num = std::to_string(abs(isoValue));
		std::string end = ".vtk";
		std::string datafile = temp + num + end;
		writer->SetFileName(datafile.c_str());
		writer->SetFileTypeToBinary();
		writer->Write();
		*/
	}

	vtkDataArray* colorArray = dataset->GetCellData()->GetArray(arrayName.c_str());
	vtkSmartPointer<vtkIdList> colorList = vtkSmartPointer<vtkIdList>::New();
	for (vtkIdType arrayId = 0; arrayId < colorArray->GetNumberOfTuples(); arrayId++)
	{
		vtkIdType currColor = static_cast<vtkIdType>(colorArray->GetTuple1(arrayId));
		colorList->InsertUniqueId(currColor);
	}

	cout << "Total Regions in this dataset: " << colorList->GetNumberOfIds() << endl;
	if (logFile.is_open())
	{
		logFile << "Total Regions in this dataset: " << colorList->GetNumberOfIds() << endl;
	}

	/*
	Json::StyledStreamWriter writer;
	std::ofstream outFile;
	outFile.open("./vortex_profile_construction/fortRegionSplitting.json");
	writer.write(outFile, tree);
	outFile.close();
	*/
	return numRegions;
}

std::vector<double> readLambda2Steps(std::string filepath, std::ofstream& logFile)
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

		cout << index << lambda2Value << endl;

		if (logFile.is_open())
		{
			logFile << index << lambda2Value << endl;
		}

		lambda2Values.push_back(lambda2Value);
	}

	return lambda2Values;
}


void createLookupTable(vtkLookupTable* lut, vtkNamedColors* colors, int numberOfRegions)
{
	lut->SetNumberOfTableValues(std::max(numberOfRegions + 1, 10));
	lut->Build();
	lut->SetTableValue(0, colors->GetColor4d("Black").GetData());
	lut->SetTableValue(1, colors->GetColor4d("MediumVioletRed").GetData());
	lut->SetTableValue(2, colors->GetColor4d("OrangeRed").GetData());
	lut->SetTableValue(3, colors->GetColor4d("DarkKhaki").GetData());
	lut->SetTableValue(4, colors->GetColor4d("Indigo").GetData());
	lut->SetTableValue(5, colors->GetColor4d("DarkGreen").GetData());
	lut->SetTableValue(6, colors->GetColor4d("DarkBlue").GetData());
	lut->SetTableValue(7, colors->GetColor4d("SaddleBrown").GetData());
	lut->SetTableValue(8, colors->GetColor4d("AntiqueWhite").GetData());
	lut->SetTableValue(9, colors->GetColor4d("Gray").GetData());

	// If the number of regions is larger than the number of specified colors,
	// generate some random colors.
	// Note: If a Python version is written, it is probably best to use
	//       vtkMinimalStandardRandomSequence in it and here, to ensure
	//       that the random number generation is the same.

	if (numberOfRegions > 9)
	{
		std::mt19937 mt(4355412); // Standard mersenne_twister_engine
		std::uniform_real_distribution<double> distribution(.1, 0.5);
		for (auto i = 10; i < numberOfRegions; ++i)
		{
			lut->SetTableValue(i, distribution(mt), distribution(mt),
				distribution(mt), 1.0);
		}
	}

	lut->SetTableValue(numberOfRegions, colors->GetColor4d("DarkRed").GetData());
}


int main(int argc, char* argv[])
{
	vtkSMPTools::SetBackend("STDThread");
	// vtkSMPTools::Initialize(8);
	int num_threads = vtkSMPTools::GetEstimatedNumberOfThreads();
	cout << "Total threads: " << num_threads << endl;
	float vsf;
	if (argc == 6)
	{
		vsf = std::stof(argv[5]);
	}
	else
	{
		vsf = 3.5;
	}

	if (argc < 5)
	{
		std::cerr << "Please specify the input filename, csv filename, output filename and json filename." << endl;
		exit(1);
	}

	vtkNew<vtkNamedColors> colors;

	// std::string datafile = "./initial_value_problem/fortRegionGrowing_float_bin.vtk";
	std::string datafile = argv[1];
	vtkNew<vtkDataSetReader> dataReader;
	dataReader->SetFileName(datafile.c_str());
	dataReader->Update();
	dataReader->GetOutput()->GetPointData()->SetActiveScalars("lambda2-vtk");

	vtkSmartPointer<vtkDataArray> arr = dataReader->GetOutput()->GetCellData()->GetArray("GlobalCellIds");
	vtkSmartPointer<vtkDataSet> dataset;

	if (arr != nullptr)
	{
		dataReader->GetOutput()->GetCellData()->RemoveArray("GlobalCellIds");
		dataReader->GetOutput()->GetPointData()->RemoveArray("GlobalPointIds");
		dataReader->GetOutput()->GetPointData()->RemoveArray("vtkGhostType");
		dataReader->Update();
	}

	// Generate new cell Ids for the candidate regions
	vtkNew<vtkGenerateGlobalIds> generator;
	generator->SetInputData(dataReader->GetOutput());
	generator->Update();
	dataset = vtkDataSet::SafeDownCast(generator->GetOutput());

	dataset->GetCellData()->SetActiveGlobalIds("");
	dataset->GetPointData()->SetActiveGlobalIds("");

	std::ofstream logFile;
	logFile.open("fortLogs.txt");
	// std::string stepsFile = "./initial_value_problem/fortFull_bin.csv";
	std::string stepsFile = argv[2];
	std::vector<double> lambda2Values = readLambda2Steps(stepsFile, logFile);

	/*
	// Extract the largest region
	vtkNew<vtkConnectivityFilter> extractLargestRegion;
	extractLargestRegion->SetInputData(dataReader->GetOutput());
	extractLargestRegion->SetExtractionModeToLargestRegion();
	extractLargestRegion->Update();
	vtkSmartPointer<vtkDataSet> bigDataset = vtkDataSet::SafeDownCast(extractLargestRegion->GetOutput());

	// Generate new cell Ids for the candidate regions
	vtkNew<vtkGenerateGlobalIds> generator;
	generator->SetInputData(bigDataset);
	generator->Update();
	vtkSmartPointer<vtkDataSet> dataset = vtkDataSet::SafeDownCast(generator->GetOutput());
	// cout << *dataset << endl;
	*/

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	Json::Value tree;
	int numRegions = colorRegions(dataset, lambda2Values, vsf, num_threads, logFile, tree);
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
	if (logFile.is_open())
	{
		logFile << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
		logFile.close();
	}

	vtkNew<vtkFloatArray> imageColor;
	imageColor->SetNumberOfComponents(1);
	imageColor->SetNumberOfTuples(dataset->GetNumberOfPoints());
	imageColor->SetName("colors");
	// Set the color for all points to 1
	for (vtkIdType ptIdx = 0; ptIdx < dataset->GetNumberOfPoints(); ptIdx++) {
		double nullId[1] = { 1 };
		imageColor->SetTuple(ptIdx, nullId);
	}
	dataset->GetPointData()->AddArray(imageColor);

	// Clean-up
	if (dataset->GetPointData()->HasArray("RegionId"))
	{
		dataset->GetPointData()->RemoveArray("RegionId");
	}
	if (dataset->GetPointData()->HasArray("vtkGhostType"))
	{
		dataset->GetPointData()->RemoveArray("vtkGhostType");
	}
	if (dataset->GetCellData()->HasArray("RegionId"))
	{
		dataset->GetCellData()->RemoveArray("RegionId");
	}

	vtkNew<vtkDataSetWriter> writer;
	writer->SetInputData(dataset);
	// std::string temp = "./vortex_profile_construction/fortRegionSplitting.vtk";
	std::string temp = argv[3];
	writer->SetFileName(temp.c_str());
	writer->SetFileTypeToBinary();
	writer->Write();

	std::string jsonFileName = argv[4];
	Json::StyledStreamWriter jsonWriter;
	std::ofstream outFile;
	outFile.open(jsonFileName);
	jsonWriter.write(outFile, tree);
	outFile.close();

	vtkNew<vtkLookupTable> lut;
	createLookupTable(lut, colors, numRegions);

	vtkNew<vtkDataSetMapper> regionMapper;
	// regionMapper->ScalarVisibilityOff();
	regionMapper->SetInputData(dataset);
	regionMapper->SetScalarModeToUseCellData();
	regionMapper->SetScalarRange(0, numRegions - 1);
	regionMapper->SetLookupTable(lut);
	regionMapper->Update();
	vtkNew<vtkActor> regionActor;
	// regionActor->GetProperty()->SetColor(colors->GetColor3d("Chartreuse").GetData());
	regionActor->GetProperty()->SetOpacity(0.5);
	regionActor->SetMapper(regionMapper);

	/*
	vtkNew<vtkDataSetMapper> cMapper;
	// cMapper->ScalarVisibilityOff();
	cMapper->SetInputData(disconnectedIsosurfaces->GetOutput());
	cMapper->SetScalarRange(0, numRegions - 1);
	cMapper->SetLookupTable(lut);
	cMapper->Update();
	vtkNew<vtkActor> cActor;
	// cActor->GetProperty()->SetColor(colors->GetColor3d("SkyBlue").GetData());
	// cActor->GetProperty()->SetOpacity(0.5);
	cActor->SetMapper(cMapper);
	*/

	vtkNew<vtkOutlineFilter> outline;
	outline->SetInputConnection(dataReader->GetOutputPort());
	outline->Update();
	vtkNew<vtkDataSetMapper> outlineMapper;
	outlineMapper->SetInputConnection(outline->GetOutputPort());
	vtkNew<vtkActor> outlineActor;
	outlineActor->SetMapper(outlineMapper);
	outlineActor->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());
	outlineActor->GetProperty()->SetLineWidth(3.0);

	vtkColor3d axesActorColor = colors->GetColor3d("Black");
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
	ren->AddActor(regionActor);
	// ren->AddActor(cActor);
	ren->SetBackground(colors->GetColor3d("White").GetData());
	ren->GetActiveCamera()->SetPosition(0, -1, 0);
	ren->GetActiveCamera()->SetViewUp(0, 0, 1);
	ren->ResetCamera();
	// ren->GetActiveCamera()->SetPosition(0, 0, 1);
	ren->GetActiveCamera()->Azimuth(-30);
	ren->GetActiveCamera()->Elevation(15);
	ren->GetActiveCamera()->Zoom(1.65);

	vtkNew<vtkRenderWindow> ren_win;
	ren_win->SetSize(1200, 600);
	ren_win->SetWindowName("Vortex Core w. Q, delta and lambda Criteria");
	ren_win->AddRenderer(ren);

	vtkNew<vtkRenderWindowInteractor> iren;
	vtkNew<vtkInteractorStyleTrackballCamera> irenStyle;
	iren->SetInteractorStyle(irenStyle);
	iren->SetRenderWindow(ren_win);

	ren_win->Render();
	iren->Start();

	return EXIT_SUCCESS;
}