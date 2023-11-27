// In this version, we use both geometric and physical criteria to simplify regions
// In this version, we remove the unnecessary cell data after processing the regions.
// In this version, we add additional geometric criteria in the vortex profile.

#include <vtkActor.h>
#include <vtkAppendFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkCamera.h>
#include <vtkCleanPolyData.h>
#include <vtkCharArray.h>
#include <vtkCellData.h>
#include <vtkColorSeries.h>
#include <vtkCompositeDataWriter.h>
#include <vtkConnectivityFilter.h>
#include <vtkContourFilter.h>
#include <vtkCubeAxesActor.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>
#include <vtkGaussianKernel.h>
#include <vtkGenerateGlobalIds.h>
#include <vtkGeometryFilter.h>
#include <vtkIntArray.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkjsoncpp/json/json.h>
#include <vtkLookupTable.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkOBBTree.h>
#include <vtkOutlineFilter.h>
#include <vtkPointData.h>
#include <vtkPointDataToCellData.h>
#include <vtkPointInterpolator.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyLine.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSMPTools.h>
#include <vtkStreamTracer.h>
#include <vtkStaticPointLocator.h>
#include <vtkSTLWriter.h>
#include <vtkStripper.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkTextProperty.h>
#include <vtkThreshold.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLMultiBlockDataWriter.h>

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
			double pt[3] = { round_up(results.at(0), 3), round_up(results.at(1), 3),
				round_up(results.at(2), 3) };
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

void extractSurface(vtkDataSet* dataset, vtkPolyData* surface,
	vtkPolyData* skeleton, int& Level, int& RegionId)
{
	vtkNew<vtkNamedColors> colors;
	cout << "Extracting Geometry...";
	vtkSmartPointer<vtkGeometryFilter> surfaceFilter = vtkSmartPointer<vtkGeometryFilter>::New();
	surfaceFilter->SetInputData(dataset);
	// surfaceFilter->Update();
	cout << " Done" << endl;
	cout << "Start surface smoothing... ";
	vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smoothFilter->SetInputConnection(surfaceFilter->GetOutputPort());
	smoothFilter->SetNumberOfIterations(10);
	smoothFilter->SetRelaxationFactor(0.25);
	smoothFilter->FeatureEdgeSmoothingOff();
	smoothFilter->BoundarySmoothingOn();
	smoothFilter->Update();
	cout << " Done" << endl;
	vtkPolyData* surfaceData = smoothFilter->GetOutput();

	if (surfaceData->GetNumberOfPoints() == 0 || surfaceData->GetNumberOfCells() == 0)
	{
		surface = nullptr;
		skeleton = nullptr;
		return;
	}

	std::string arr_name = "RegionIds_" + std::to_string(Level);
	std::vector<std::string> arrsToDelte;
	for (int i = 0; i < surfaceData->GetCellData()->GetNumberOfArrays(); i++)
	{
		if (surfaceData->GetCellData()->GetArrayName(i) != arr_name)
		{
			arrsToDelte.push_back(surfaceData->GetCellData()->GetArrayName(i));
		}
	}
	for (int i = 0; i < arrsToDelte.size(); i++)
	{
		surfaceData->GetCellData()->RemoveArray(arrsToDelte.at(i).c_str());
	}

	surface->DeepCopy(surfaceData);

	std::string filename = "mesh";
	std::string ext = ".stl";
	std::string mesh_filename_wExt = filename + ext;
	// First write the extracted surface to an STL file
	std::string mesh_filepath = mesh_filename_wExt;
	cout << "Writing surface data...";
	vtkNew<vtkSTLWriter> stlWriter;
	stlWriter->SetFileName(mesh_filepath.c_str());
	stlWriter->SetInputData(surfaceData);
	stlWriter->SetFileTypeToBinary();
	stlWriter->Write();
	cout << " Done" << endl;
	cout << "Extracting skeleton...";
	std::string skel_exe_path = "Extract_Skeleton.exe";
	std::string skel_output_path = filename + ".txt";
	std::string command = skel_exe_path + " \"" + mesh_filepath + "\" \"" + skel_output_path + "\"";
	system(command.c_str());
	cout << " Done" << endl;
	cout << "Reading skel info...";
	vtkSmartPointer<vtkPolyData> skel = vtkSmartPointer<vtkPolyData>::New();
	readSkelInfo(skel_output_path, skel);
	cout << " Done" << endl;
	cout << "Stripping skeleton....";
	if (skel->GetNumberOfCells() == 0 ||
		skel->GetNumberOfPoints() == 0)
	{
		surface = nullptr;
		skeleton = nullptr;
		return;
	}

	vtkNew<vtkCleanPolyData> clean;
	clean->SetInputData(skel);
	clean->Update();
	// cout << *clean->GetOutput() << endl;

	vtkNew<vtkStripper> stripper;
	stripper->SetInputData(clean->GetOutput());
	stripper->JoinContiguousSegmentsOn();
	stripper->Update();
	cout << " Done" << endl;
	cout << "Smoothing skeleton...";
	smoothFilter = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
	smoothFilter->SetInputData(stripper->GetOutput());
	smoothFilter->SetNumberOfIterations(10);
	smoothFilter->SetRelaxationFactor(0.5);
	smoothFilter->FeatureEdgeSmoothingOff();
	smoothFilter->BoundarySmoothingOn();
	smoothFilter->Update();
	cout << " Done" << endl;
	cout << "Populating data on skeleton...";
	vtkSmartPointer<vtkGaussianKernel> kernel = vtkSmartPointer<vtkGaussianKernel>::New();
	kernel->SetKernelFootprintToNClosest();
	kernel->SetNumberOfPoints(8);
	
	vtkSmartPointer<vtkPointInterpolator> interpolator = vtkSmartPointer<vtkPointInterpolator>::New();
	interpolator->SetInputData(smoothFilter->GetOutput());
	interpolator->SetSourceData(dataset);
	interpolator->SetKernel(kernel);
	interpolator->SetNullPointsStrategyToClosestPoint();
	interpolator->Update();
	skeleton->DeepCopy(interpolator->GetPolyDataOutput());
	cout << " Done" << endl;
	// Add regionIds array
	cout << "Assigning Ids to skeleton...";
	vtkNew<vtkIdTypeArray> skelRegionIdsArray;
	skelRegionIdsArray->SetName(arr_name.c_str());
	skelRegionIdsArray->SetNumberOfTuples(skeleton->GetNumberOfCells());
	for (int i = 0; i < skeleton->GetNumberOfCells(); i++)
	{
		skelRegionIdsArray->SetTuple1(i, RegionId);
	}
	skeleton->GetCellData()->AddArray(skelRegionIdsArray);
	cout << " Done" << endl;
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
	regMapper->SetInputData(dataset);
	regMapper->ScalarVisibilityOff();
	regMapper->Update();
	vtkNew<vtkActor> regActor;
	regActor->GetProperty()->SetColor(colors->GetColor3d("Blue").GetData());
	regActor->SetMapper(regMapper);
	// regActor->GetProperty()->SetLineWidth(2.0);
	regActor->GetProperty()->SetOpacity(0.1);
	ren1->AddActor(regActor);

	vtkNew<vtkDataSetMapper> smallRegMapper;
	smallRegMapper->SetInputData(surface);
	smallRegMapper->ScalarVisibilityOff();
	smallRegMapper->Update();
	vtkNew<vtkActor> smallRegActor;
	smallRegActor->GetProperty()->SetColor(colors->GetColor3d("Green").GetData());
	smallRegActor->GetProperty()->SetLineWidth(2.0);
	smallRegActor->GetProperty()->SetOpacity(0.25);
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
	return;
}

std::vector<int> linspace(int start_in, int end_in, int num_in)
{
	std::vector<int> linspaced;

	if (end_in < num_in)
	{
		num_in = end_in;
	}

	if (num_in == 0)
	{
		return linspaced;
	}

	linspaced.push_back(start_in);

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

double CalculateGaussCurvature(vtkDataSet* polyData)
{
	// cout << "Curveature Start..." << endl;
	int numCells = polyData->GetNumberOfCells();
	double sumAngles = 0.0;
	int totalSteps = 0;
	for (int cellId = 0; cellId < numCells; cellId += 1)
	{
		vtkCell* polyLine = polyData->GetCell(cellId);
		vtkPoints* points = polyLine->GetPoints();
		int numPoints = points->GetNumberOfPoints();
		totalSteps += numPoints;
		// cout << "NumPts: " << numPoints << endl;
		if (numPoints == 1)
		{
			continue;
		}
		std::vector<int> steps = linspace(0, numPoints, numPoints);
		// cout << "NumStps: " << steps.size() << endl;
		if (steps.size() < 3 || (numCells == 1 && numPoints == 1))
		{
			return 0.;
		}
		// totalSteps += steps.size();
		// cout << totalSteps << " Steps: " << numPoints << endl;
		// Compute the tangent vectors
		// std::vector<double> tangentX(steps.size()), tangentY(steps.size()), tangentZ(steps.size());
		// for (int i = 1; i < steps.size() - 1; i++) {

		for (int i = 1; i < numPoints - 1; i += 2) {
			double prevPt[3], currPt[3], nextPt[3];
			points->GetPoint(i - 1, prevPt);
			points->GetPoint(i, currPt);
			points->GetPoint(i + 1, nextPt);
			// cout << currPt[0] << " " << currPt[1] << " " << currPt[2];
			double dX1 = currPt[0] - prevPt[0];
			double dY1 = currPt[1] - prevPt[1];
			double dZ1 = currPt[2] - prevPt[2];
			double dX2 = nextPt[0] - currPt[0];
			double dY2 = nextPt[1] - currPt[1];
			double dZ2 = nextPt[2] - currPt[2];
			// cout << steps[i - 1] << " " << steps[i] << " " << steps[i + 1] << endl;
			// tangentX[i] = dY1 * dZ2 - dZ1 * dY2;
			// tangentY[i] = dZ1 * dX2 - dX1 * dZ2;
			// tangentZ[i] = dX1 * dY2 - dY1 * dX2;
			double tangent1[3] = { dX1, dY1, dZ1 };
			double tangent2[3] = { dX2, dY2, dZ2 };
			vtkMath::Normalize(tangent1);
			vtkMath::Normalize(tangent2);
			// double mag = sqrt(tangent[0] * tangent[0] + tangent[1] * tangent[1] + tangent[2] * tangent[2]);
			// double mag = sqrt(tangentX[i] * tangentX[i] + tangentY[i] * tangentY[i] + tangentZ[i] * tangentZ[i]);
			// tangentX[i] /= mag;
			// tangentY[i] /= mag;
			// tangentZ[i] /= mag;
			double dotProd = tangent1[0] * tangent2[0] + tangent1[1] * tangent2[1] + tangent1[2] * tangent2[2];
			double angle = acos(dotProd);

			// cout << tangent1[0] << " " << tangent1[1] << " " << tangent1[2] << endl;
			// cout << tangent2[0] << " " << tangent2[1] << " " << tangent2[2] << endl;

			if (i == 1)
			{
				angle = 0;
			}

			if (std::isnan(angle))
			{
				continue;
			}
			sumAngles += angle;

			// cout << " dotProd: " << dotProd << " sumAngles: " << sumAngles << " myNumPoints: " << myNumPoints << endl;
		}
	}
	/*
	// Compute the curvature using the tangent vectors

	for (int i = 1; i < numPoints - 1; i++) {
		double dotProd = tangentX[i - 1] * tangentX[i] + tangentY[i - 1] * tangentY[i] + tangentZ[i - 1] * tangentZ[i];
		cout << dotProd;

		cout << " sumAngles: " << sumAngles << endl;
	}
	*/
	// Compute the Gauss curvature
	double gaussCurvature = sumAngles / (pow(totalSteps, 1)) * 100;
	// cout << "Curveature End..." << endl;
	return gaussCurvature;
}

std::vector<double> inspectSkeletonOrientation(vtkPolyData* skeleton)
{
	std::vector<double> geometry_info;

	int numPts, numCells = skeleton->GetNumberOfCells(), numVec = 0, totalPts = 0;
	double i_unit[3] = { 1.0, 0.0, 0.0 };
	double j_unit[3] = { 0.0, 1.0, 0.0 };
	double k_unit[3] = { 0.0, 0.0, 1.0 };
	double prevPoint[3], currPoint[3], currVec[3];
	double length, lineLength, xlength, yLength, zLength;
	double totallength = 0.0, totalxLength = 0.0, totalyLength = 0.0, totalzLength = 0.0;
	double xDir, yDir, zDir, maxDir;
	vtkIdType pointId;

	for (int i = 0; i < numCells; i++)
	{
		lineLength = 0.0;
		xlength = 0.0;
		yLength = 0.0;
		zLength = 0.0;
		vtkCell* polyLine = skeleton->GetCell(i);
		vtkPoints* points = polyLine->GetPoints();
		vtkIdList* pointIds = polyLine->GetPointIds();

		numPts = points->GetNumberOfPoints();
		totalPts += numPts;
		for (int j = 0; j < numPts - 1; j++)
		{
			numVec++;
			pointId = pointIds->GetId(j);
			points->GetPoint(j, prevPoint);
			points->GetPoint(j + 1, currPoint);
			length = sqrt(vtkMath::Distance2BetweenPoints(prevPoint, currPoint));
			lineLength += length;
			currVec[0] = currPoint[0] - prevPoint[0];
			currVec[1] = currPoint[1] - prevPoint[1];
			currVec[2] = currPoint[2] - prevPoint[2];
			vtkMath::Normalize(currVec);
			xDir = vtkMath::Dot(currVec, i_unit);
			xDir = std::sqrt(xDir*xDir);	// Get absolute value
			yDir = vtkMath::Dot(currVec, j_unit);
			yDir = std::sqrt(yDir*yDir);
			zDir = vtkMath::Dot(currVec, k_unit);
			zDir = std::sqrt(zDir*zDir);
			xlength += xDir;
			yLength += yDir;
			zLength += zDir;
		}

		totallength += lineLength;
		totalxLength += xlength;
		totalyLength += yLength;
		totalzLength += zLength;
	}

	double xPerc = (totalxLength / numVec);
	double yPerc = (totalyLength / numVec);
	double zPerc = (totalzLength / numVec);
	geometry_info.push_back(xPerc);
	geometry_info.push_back(yPerc);
	geometry_info.push_back(zPerc);

	double curvature = CalculateGaussCurvature(skeleton);
	double gauss_curvature = curvature * (1 - xPerc) * (yPerc) * (zPerc);
	double bbox_length = skeleton->GetLength();

	double 	cor[3], max[3], mid[3], min[3], size[3];
	vtkOBBTree::ComputeOBB(skeleton->GetPoints(), cor, max, mid, min, size);

	double p1[3] = { cor[0] + max[0], cor[1] + max[1], cor[2] + max[2] };
	double p2[3] = { cor[0] + mid[0], cor[1] + mid[1], cor[2] + mid[2] };
	double p3[3] = { cor[0] + min[0], cor[1] + min[1], cor[2] + min[2] };
	double a = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
	double b = sqrt(vtkMath::Distance2BetweenPoints(p1, p3));
	double c = sqrt(vtkMath::Distance2BetweenPoints(p2, p3));
	// Pythagorous Theorem
	bbox_length = sqrt(pow(sqrt(pow(a, 2) + pow(b, 2)), 2) + pow(c, 2));

	double skelRatio = bbox_length / totallength;
	double newSkelLen = totallength / skelRatio;
	double newSkelRatio = newSkelLen / totallength;
	double normal_length = totallength / std::pow(numVec, 0.25);
	double newCurv = (gauss_curvature * newSkelRatio * normal_length);
	skelRatio = totallength / bbox_length;
	geometry_info.push_back(curvature);
	geometry_info.push_back(newCurv);
	geometry_info.push_back(totallength);
	geometry_info.push_back(skelRatio);
	return geometry_info;
}

std::string buildVortexProfileVectors(std::vector<double>& vortexProfile, vtkDataSet* newRegDataset,
	vtkPolyData* skeleton)
{
	/*
	0. Avg(lambda2)					1. Avg(lambdaci)
	2. Avg(Q)						3. Avg(Delta)
	4. Avg(Divergence)				5. Avg(Spanwise)
	6. Size							7. Avg(Vorticity)
	8. Avg(Enstrophy)				9. Avg(Velocity)
	10. Avg(Acceleration)			11. Norm(Jacobian)
	12. Curvature					13. FullCurvature
	14. Streamwise Direction(S_t)	15. Spanwise direction(S_p)
	16. Vertical Direction(S_v)		17. Length(L)
	18. rho
	*/

	std::string name = "";
	double arraySum = 0;
	double arrayAvg;
	int nTuples = newRegDataset->GetNumberOfPoints();
	cout << "0. Avg(lambda2)" << endl;
	// 0. Avg(lambda2)
	vtkDataArray* arr = newRegDataset->GetPointData()->GetArray("lambda2-vtk");
	for (vtkIdType i = 0; i < nTuples; i++)
	{
		double value = arr->GetTuple1(i);
		arraySum += value;
	}
	arrayAvg = (arraySum / nTuples);
	if (std::isnan(arrayAvg))
	{
		arrayAvg = 0;
	}
	name += "_Lambda2=" + std::to_string(std::abs(arrayAvg));
	vortexProfile.push_back(arrayAvg);
	cout << "1. Avg(lambdaci)" << endl;
	// 1. Avg(lambdaci)
	arr = newRegDataset->GetPointData()->GetArray("lambdaci-vtk");
	arraySum = 0;
	for (vtkIdType i = 0; i < nTuples; i++)
	{
		double value = arr->GetTuple1(i);
		arraySum += value;
	}
	arrayAvg = (arraySum / nTuples);
	if (std::isnan(arrayAvg))
	{
		arrayAvg = 0;
	}
	name += "_Lambdaci=" + std::to_string(std::abs(arrayAvg));
	vortexProfile.push_back(arrayAvg);
	cout << "2. Avg(Q)" << endl;
	// 2. Avg(Q)
	arr = newRegDataset->GetPointData()->GetArray("q-vtk");
	arraySum = 0;
	for (vtkIdType i = 0; i < nTuples; i++)
	{
		double value = arr->GetTuple1(i);
		arraySum += value;
	}
	arrayAvg = (arraySum / nTuples);
	if (std::isnan(arrayAvg))
	{
		arrayAvg = 0;
	}
	name += "_Q=" + std::to_string(std::abs(arrayAvg));
	vortexProfile.push_back(arrayAvg);
	cout << "3. Avg(Delta)" << endl;
	// 3. Avg(Delta)
	arr = newRegDataset->GetPointData()->GetArray("delta-vtk");
	arraySum = 0;
	for (vtkIdType i = 0; i < nTuples; i++)
	{
		double value = arr->GetTuple1(i);
		arraySum += value;
	}
	arrayAvg = (arraySum / nTuples);
	if (std::isnan(arrayAvg))
	{
		arrayAvg = 0;
	}
	name += "_Delta=" + std::to_string(std::abs(arrayAvg));
	vortexProfile.push_back(arrayAvg);
	cout << "4. Avg(Divergence)" << endl;
	// 4. Avg(Divergence)
	arr = newRegDataset->GetPointData()->GetArray("divergence-vtk");
	arraySum = 0;
	for (vtkIdType i = 0; i < nTuples; i++)
	{
		double value = arr->GetTuple1(i);
		arraySum += value;
	}
	arrayAvg = (arraySum / nTuples);
	if (std::isnan(arrayAvg))
	{
		arrayAvg = 0;
	}
	name += "_Divergence=" + std::to_string(arrayAvg);
	vortexProfile.push_back(arrayAvg);
	cout << "5. Avg(oyf)" << endl;
	// 5. Avg(oyf)
	if (newRegDataset->GetPointData()->HasArray("oyf"))
	{
		arr = newRegDataset->GetPointData()->GetArray("oyf");
		arraySum = 0;

		for (vtkIdType i = 0; i < nTuples; i++)
		{
			double value = arr->GetTuple1(i);

			if (value < 0)
			{
				value = value * 0.1;
			}
			value = std::pow(value, 2);
			arraySum += value;
		}

		arrayAvg = std::sqrt((arraySum / nTuples));
		if (std::isnan(arrayAvg))
		{
			arrayAvg = 0;
		}
		name += "_Oyf=" + std::to_string(arrayAvg);
		vortexProfile.push_back(arrayAvg);
	}
	cout << "6. Size" << endl;
	// 6. Size
	double size = static_cast<double>(newRegDataset->GetNumberOfCells());
	name += "_Size=" + std::to_string(size);
	vortexProfile.push_back(size);
	cout << "7. Avg(Vorticity)" << endl;
	// 7. Avg(Vorticity)
	// 8. Avg(Enstrophy)
	arr = newRegDataset->GetPointData()->GetArray("vorticity-vtk");
	double magSum = 0;
	double ensSum = 0;
	for (vtkIdType i = 0; i < nTuples; i++)
	{
		double values[3];
		arr->GetTuple(i, values);
		double mag = vtkMath::Norm(values);
		magSum += mag;
		double ens = 0.5 * vtkMath::SquaredNorm(values);
		ensSum += ens;
	}
	arrayAvg = (magSum / nTuples);
	if (std::isnan(arrayAvg))
	{
		arrayAvg = 0;
	}
	name += "_Vorticity=" + std::to_string(arrayAvg);
	vortexProfile.push_back(arrayAvg);
	cout << "8. Avg(Enstrophy)" << endl;
	arrayAvg = (ensSum / nTuples);
	if (std::isnan(arrayAvg))
	{
		arrayAvg = 0;
	}
	name += "_Enstrophy=" + std::to_string(arrayAvg);
	vortexProfile.push_back(arrayAvg);
	cout << "9. Avg(Velocity)" << endl;
	// 9. Avg(Velocity)
	arr = newRegDataset->GetPointData()->GetArray("velocity");
	for (vtkIdType i = 0; i < nTuples; i++)
	{
		double values[3];
		arr->GetTuple(i, values);
		double mag = vtkMath::Norm(values);
		arraySum += mag;
	}
	arrayAvg = (arraySum / nTuples);
	if (std::isnan(arrayAvg))
	{
		arrayAvg = 0;
	}
	name += "_Velocity=" + std::to_string(arrayAvg);
	vortexProfile.push_back(arrayAvg);
	cout << "10. Avg(Acceleration)" << endl;
	// 10. Avg(Acceleration)
	arr = newRegDataset->GetPointData()->GetArray("acceleration-vtk");
	for (vtkIdType i = 0; i < nTuples; i++)
	{
		double values[3];
		arr->GetTuple(i, values);
		double mag = vtkMath::Norm(values);
		arraySum += mag;
	}
	arrayAvg = (arraySum / nTuples);
	if (std::isnan(arrayAvg))
	{
		arrayAvg = 0;
	}
	name += "_Acceleration=" + std::to_string(arrayAvg);
	vortexProfile.push_back(arrayAvg);
	cout << "11. Norm(Jacobian)" << endl;
	// 11. Norm(Jacobian)
	arr = newRegDataset->GetPointData()->GetArray("jacobian-vtk");
	for (vtkIdType i = 0; i < nTuples; i++)
	{
		double v[9];
		arr->GetTuple(i, v);
		double mag = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2] +
			v[3] * v[3] + v[4] * v[4] + v[5] * v[5] + v[6] * v[6] + v[7] * v[7] + v[8] * v[8]);
		arraySum += mag;
	}
	arrayAvg = (arraySum / nTuples);
	if (std::isnan(arrayAvg))
	{
		arrayAvg = 0;
	}
	name += "_Jacobian=" + std::to_string(arrayAvg);
	vortexProfile.push_back(arrayAvg);

	std::vector<double> geometry_info = inspectSkeletonOrientation(skeleton);
	cout << "12. Curvature" << endl;
	// 12. Curvature
	if (std::isnan(geometry_info[3]))
	{
		geometry_info[3] = 0.0;
	}
	name += "_Curvature=" + std::to_string(geometry_info[3]);
	vortexProfile.push_back(geometry_info[3]);
	cout << "13. Full Curvature" << endl;
	// 13. Full Curvature
	if (std::isnan(geometry_info[4]))
	{
		geometry_info[4] = 0.0;
	}
	name += "_FullCurvature=" + std::to_string(geometry_info[4]);
	vortexProfile.push_back(geometry_info[4]);
	cout << "14. Streamwise Direction(S_t)" << endl;
	// 14. Streamwise Direction(S_t)
	if (std::isnan(geometry_info[0]))
	{
		geometry_info[0] = 0.0;
	}
	name += "_St=" + std::to_string(geometry_info[0]);
	vortexProfile.push_back(geometry_info[0]);
	cout << "15. Spanwise direction(S_p)" << endl;
	// 15. Spanwise direction(S_p)
	if (std::isnan(geometry_info[1]))
	{
		geometry_info[1] = 0.0;
	}
	name += "_Sp=" + std::to_string(geometry_info[1]);
	vortexProfile.push_back(geometry_info[1]);
	cout << "16. Vertical Direction(S_v)" << endl;
	// 16. Vertical Direction(S_v)
	if (std::isnan(geometry_info[2]))
	{
		geometry_info[2] = 0.0;
	}
	name += "_Sv=" + std::to_string(geometry_info[2]);
	vortexProfile.push_back(geometry_info[2]);
	cout << "17. Length(L)" << endl;
	// 17. Length(L)
	if (std::isnan(geometry_info[5]))
	{
		geometry_info[5] = 0.0;
	}
	name += "_Length=" + std::to_string(geometry_info[5]);
	vortexProfile.push_back(geometry_info[5]);
	cout << "18. rho" << endl;
	// 18. rho
	if (std::isnan(geometry_info[6]))
	{
		geometry_info[6] = 0.0;
	}
	name += "_Rho=" + std::to_string(geometry_info[6]);
	vortexProfile.push_back(geometry_info[6]);

	return name;
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

int regionSimplification(vtkDataSet* regDataset, vtkDataSet* outputGrid, std::string critArrName)
{
	vtkSmartPointer<vtkPointDataToCellData> ptToCell = vtkSmartPointer<vtkPointDataToCellData>::New();
	ptToCell->PassPointDataOn();
	ptToCell->SetInputData(regDataset);
	ptToCell->Update();
	vtkSmartPointer<vtkUnstructuredGrid> cellData = vtkUnstructuredGrid::SafeDownCast(ptToCell->GetOutput());
	/*
	vtkSmartPointer<vtkUnstructuredGrid> cellData = vtkSmartPointer<vtkUnstructuredGrid>::New();
	cellData->DeepCopy(ptToCell->GetUnstructuredGridOutput());
	*/

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
	return *numRegs;
}


void updateTreeJson(vtkDataSet* dataset, Json::Value& tree, vtkPolyData* surface, vtkPolyData* skeleton,
	int& Level, int& RegionId)
{
	std::string regionName = std::to_string(RegionId);
	Json::Value &nodes = tree["children"];

	for (unsigned i = 0; i < nodes.size(); i++)
	{
		Json::Value &node = nodes[i];
		std::string nodeName = node["name"].asString();
		std::size_t pos1 = nodeName.find("_");
		std::string regionId = nodeName.substr(0, pos1);
		std::size_t pos2 = nodeName.find("_", pos1 + 1);
		std::string level = nodeName.substr(pos1 + 1, pos2 - pos1 - 1);
		std::size_t pos3 = nodeName.find("_", pos2 + 1);
		std::string parentId = nodeName.substr(pos2 + 1, pos3 - pos2 - 1);
		std::size_t pos4 = nodeName.find("_", pos3 + 1);
		std::string lambda2 = nodeName.substr(pos3 + 1, pos4 - pos3 - 1);
		// cout << nodeName << " " << regionId << " " << level << " " << parentId << " " << lambda2 << endl;
		if (regionName == regionId && level == std::to_string(Level))
		{
			Json::Value newNode;
			std::string newName = regionId + "_" + level + "_" + parentId + "_Isovalue=" + lambda2;
			cout << newName << endl;
			vtkSmartPointer<vtkUnstructuredGrid> simple_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
			std::string criteriaArrName = "lambda2-vtk";
			cout << "Before simplification... " << endl;
			int numRegs = regionSimplification(dataset, simple_grid, criteriaArrName);
			cout << "Simplification done... " << endl;
			if (simple_grid->GetNumberOfCells() == 0 ||
				simple_grid->GetNumberOfPoints() == 0)
			{
				continue;
			}

			cout << "Start extract surface... " << endl;
			extractSurface(simple_grid, surface, skeleton, Level, RegionId);

			if (skeleton == nullptr || surface == nullptr)
			{
				continue;
			}

			cout << "Surface extraction done... " << endl;
			std::vector<double> vortexProfile;
			newName += buildVortexProfileVectors(vortexProfile, simple_grid, skeleton);

			/*
			int nCriteria = vortexProfile.size();
			// Add the vortex profile to the tree node name for later retrieval
			for (int ii = 0; ii < nCriteria; ii++)
			{
				newName += "_" + std::to_string(vortexProfile[ii]);
			}
			*/

			node["name"] = newName;
			return;
		}
		else
		{
			updateTreeJson(dataset, nodes[i], surface, skeleton, Level, RegionId);
		}
	}
	/*
	if (*numRegs <= 1)
	{
		return *numRegs;
	}

	vtkSmartPointer<vtkPolyData> psource = vtkSmartPointer<vtkPolyData>::New();
	getSeeds(psource, regDataset, 5);

	vtkSmartPointer<vtkStreamTracer> SlTracer = vtkSmartPointer<vtkStreamTracer>::New();
	SlTracer->SetInputData(regDataset);
	SlTracer->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "vorticity-vtk");
	SlTracer->SetSourceData(psource);
	SlTracer->SetMaximumPropagation(100);
	SlTracer->SetIntegratorTypeToRungeKutta45();
	SlTracer->SetIntegrationDirectionToBoth();
	SlTracer->Update();

	vtkSmartPointer<vtkPolyData> SlData = vtkSmartPointer<vtkPolyData>::New();
	double lenTh = (cellData->GetLength() * 1) / 100;
	filterSmallLines(SlTracer->GetOutput(), SlData, lenTh);

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

}

std::vector<Region> divide_regions(vtkDataSet* dataset, 
	vtkIdList* newIds, vtkDataArray* tempRegColorIds, int numThreads)
{
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
	return datasets;
}

int main(int argc, char* argv[])
{
	vtkSMPTools::SetBackend("STDThread");
	// vtkSMPTools::Initialize(8);
	int num_threads = vtkSMPTools::GetEstimatedNumberOfThreads();
	cout << "Total threads: " << num_threads << endl;

	vtkNew<vtkNamedColors> colors;
	std::string tree_file_name = argv[1];
	std::ifstream tree_file(tree_file_name);
	Json::Value tree;
	tree_file >> tree;

	std::string datafile = argv[2];
	vtkNew<vtkDataSetReader> dataReader;
	dataReader->SetFileName(datafile.c_str());
	dataReader->Update();

	vtkSmartPointer<vtkDataSet> dataset = dataReader->GetOutput();
	// Find max level in dataset
	std::vector<int> levels;
	for (int i = 1; i < dataset->GetCellData()->GetNumberOfArrays(); i++)
	{
		std::string arrayName = "RegionIds_" + std::to_string(i);
		if (dataset->GetCellData()->HasArray(arrayName.c_str()))
		{
			levels.push_back(i);
		}
	}

	vtkNew<vtkMultiBlockDataSet> skeleton_data;
	vtkNew<vtkMultiBlockDataSet> surface_data;
	std::string criteriaArrName = "lambda2-vtk";
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
	int l_idx = levels.size();
	int last_lvl = levels.at(l_idx - 1);
	for (int level : levels)
	{
		/*
		if (level != 2)
		{
			continue;
		}
		*/

		std::string arr_name = "RegionIds_" + std::to_string(level);
		vtkSmartPointer<vtkDataArray> region_ids_array = dataset->GetCellData()->GetArray(arr_name.c_str());
		vtkNew<vtkIdList> unique_region_ids;
		vtkIdType numCells = region_ids_array->GetNumberOfTuples();
		for (int i = 0; i < numCells; i++)
		{
			vtkIdType id = static_cast<vtkIdType>(region_ids_array->GetTuple1(i));
			unique_region_ids->InsertUniqueId(id);
		}
		unique_region_ids->Sort();

		std::vector<Region> datasets = divide_regions(dataset,
			unique_region_ids, region_ids_array, num_threads);

		vtkSmartPointer<vtkPolyData> skeleton_PD = nullptr;
		vtkSmartPointer<vtkPolyData> surface_PD = nullptr;
		for (int datasetId = 0; datasetId < datasets.size(); datasetId++)
		{
			vtkSmartPointer<vtkPolyData> outSkelData = nullptr;
			vtkSmartPointer<vtkPolyData> outSurfData = nullptr;
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

				cout << "Thresholding " << id << " " << level << endl;
				threshold->SetLowerThreshold(id);
				threshold->SetUpperThreshold(id);
				threshold->Update();

				vtkSmartPointer<vtkDataSet> region_data = threshold->GetOutput();
				vtkSmartPointer<vtkPolyData> skeleton = vtkSmartPointer<vtkPolyData>::New();
				vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();
				// vtkSmartPointer<vtkUnstructuredGrid> simple_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
				std::string criteriaArrName = "lambda2-vtk";
				// int numRegs = regionSimplification(region_data, simple_grid, criteriaArrName);

				// if (simple_grid->GetNumberOfCells() == 0)
				// {
				//	cout << "Wrong Regions: " << numRegs << " " << level << " " << unique_region_ids[i] << " " << endl;
				// }

				// extractSurface(simple_grid, surface, skeleton, level, unique_region_ids[i]);
				updateTreeJson(region_data, tree, surface, skeleton, level, id);

				if (surface->GetNumberOfPoints() == 0)
				{
					continue;
				}
				// cout << *surface << endl;
				if (outSkelData == nullptr)
				{
					outSkelData = vtkSmartPointer<vtkPolyData>::New();
					outSkelData->DeepCopy(skeleton);
					outSurfData = vtkSmartPointer<vtkPolyData>::New();
					outSurfData->DeepCopy(surface);
				}
				else
				{
					vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
					appendFilter->AddInputData(outSkelData);
					appendFilter->AddInputData(skeleton);
					appendFilter->Update();
					outSkelData->DeepCopy(appendFilter->GetOutput());

					appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
					appendFilter->AddInputData(outSurfData);
					appendFilter->AddInputData(surface);
					appendFilter->Update();
					outSurfData->DeepCopy(appendFilter->GetOutput());
				}
			}

			if (outSkelData == nullptr || outSurfData == nullptr)
			{
				continue;
			}
			if (skeleton_PD == nullptr)
			{
				skeleton_PD = vtkSmartPointer<vtkPolyData>::New();
				skeleton_PD->DeepCopy(outSkelData);
				surface_PD = vtkSmartPointer<vtkPolyData>::New();
				surface_PD->DeepCopy(outSurfData);
			}
			else
			{
				vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
				appendFilter->AddInputData(skeleton_PD);
				appendFilter->AddInputData(outSkelData);
				appendFilter->Update();
				skeleton_PD->DeepCopy(appendFilter->GetOutput());

				appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
				appendFilter->AddInputData(surface_PD);
				appendFilter->AddInputData(outSurfData);
				appendFilter->Update();
				surface_PD->DeepCopy(appendFilter->GetOutput());
			}
		}

		skeleton_data->SetBlock(level, skeleton_PD);
		surface_data->SetBlock(level, surface_PD);
	}
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
	// std::string temp = "./surf_skel_extraction/fortXMLSkeleton.vtm";

	/*
	vtkNew<vtkXMLMultiBlockDataWriter> xmlWriter;
	xmlWriter->SetInputData(skeleton_data);
	xmlWriter->SetFileName(temp.c_str());
	xmlWriter->SetDataModeToBinary();
	xmlWriter->Write();

	xmlWriter->SetInputData(surface_data);
	temp = "./surf_skel_extraction/fortXMLSurface.vtm";
	xmlWriter->SetFileName(temp.c_str());
	xmlWriter->SetDataModeToBinary();
	xmlWriter->Write();
	*/

	vtkNew<vtkCompositeDataWriter> writer;
	writer->SetInputData(skeleton_data);
	// temp = "./surf_skel_extraction/fortSkeleton.vtk";
	std::string temp = argv[3];
	writer->SetFileName(temp.c_str());
	writer->SetFileTypeToBinary();
	writer->Write();

	writer->SetInputData(surface_data);
	temp = argv[4];
	writer->SetFileName(temp.c_str());
	writer->SetFileTypeToBinary();
	writer->Write();

	Json::StyledStreamWriter jsonwriter;
	std::ofstream outFile;
	temp = argv[5];
	outFile.open(temp);
	jsonwriter.write(outFile, tree);
	outFile.close();

	return EXIT_SUCCESS;
}