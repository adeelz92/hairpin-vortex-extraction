// In this version, we compute the dot product and get maximum direction.
// In this version, we calculate the average dot product orientation and filter based on category.
// In this version, we perfrom the oyf analysis.
// In this version, we filter based on the curvature.
// In this version, we filter based on the combination of curvature and length.
// In this version, we save skeleton alongwith corresponding surface for visualization.
// In this version, we fix the problem of curvature being different everytime.

/*
Vortex_Statistics.exe
"Path\To\Dataset_Skeleton.vtk"
"Path\To\Dataset_Surface.vtk"
0
"Path\To\Dataset_VortSkeleton.vtk"
"Path\To\Dataset_VortSurface.vtk"
*/

#include <vtkActor.h>
#include <vtkAppendFilter.h>
#include <vtkCamera.h>
#include <vtkCellData.h>
#include <vtkCleanPolyData.h>
#include <vtkCubeAxesActor.h>
#include <vtkDataSet.h>
#include <vtkDataSetMapper.h>
#include <vtkDataSetReader.h>
#include <vtkDataSetWriter.h>
#include <vtkGeometryFilter.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkOBBTree.h>
#include <vtkOutlineFilter.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSMPTools.h>
#include <vtkStripper.h>
#include <vtkTextProperty.h>
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
	linspaced.push_back(end_in - 1);

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


double getArrAvg(vtkDataArray* array)
{
	int numTuples = array->GetNumberOfTuples();
	double value = 0;
	for (int i = 0; i < numTuples; i++)
	{
		value += array->GetTuple1(i);
	}
	return (value / numTuples);
}

double normalize(double x, double min, double max) {
	return (x - min) / (max - min);
}

double normalize_reverse(double x, const double & minVal, const double & maxVal) {
	return (2.0*(minVal - x) / (maxVal - minVal)) + 1.0;
}

struct VortexInfo
{
private:
	std::string Type;
	int Count;
public:
	VortexInfo(std::string type, int count)
	{
		this->Type = type;
		this->Count = count;
	}

	void increaseCount()
	{
		this->Count++;
	}

	std::string getType()
	{
		return this->Type;
	}

	int getCount()
	{
		return this->Count;
	}
};

// Check if type is found in the vector, 
// if Yes increment the count else add a new type and return.
void incrementType(std::vector<VortexInfo>& vortexCounts, std::string & type)
{
	bool found = false;
	for (VortexInfo& vortexCount : vortexCounts)
	{
		if (vortexCount.getType().compare(type) == 0)
		{
			vortexCount.increaseCount();
			// cout << vortexCount.getType() << " " << vortexCount.getCount() << endl;
			found = true;
			break;
		}
	}

	if (!found)
	{
		// Create a new vortex type and add in the vector
		VortexInfo info(type, 1);
		vortexCounts.push_back(info);
	}
}

int inspectSkeletonOrientation(vtkPolyData* skeleton, const double & zMin, const double & zMax,
	const double & len, std::vector<VortexInfo>& vortexCounts)
{
	int numPts, numCells = skeleton->GetNumberOfCells(), numVec = 0, totalPts = 0;
	/**/
	if (numCells > 6)
	{
		cout << "Undefined" << endl;
		std::string type = "Undefined";
		incrementType(vortexCounts, type);
		return -1;
	}

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
			// cout << pointId << " " << currPoint[0] << " " << currPoint[1] << " " << currPoint[2] << endl;
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

	/**/
	vtkDataArray* oyf = skeleton->GetPointData()->GetArray("oyf-vtk");
	double oyf_avg = 0.0;
	double numTpls = oyf->GetNumberOfTuples();
	for (int i = 0; i < numTpls; i++)
	{
		double oyf_value = oyf->GetTuple1(i);
		double nor = normalize_reverse(skeleton->GetPoint(i)[2], zMin, zMax);
		/*
		if (oyf_value < 0)
		{
			oyf_value = 0.1 * oyf_value;
		}
		oyf_value = pow(oyf_value, 2);
		*/
		oyf_avg += oyf_value;// *nor;
	}
	oyf_avg = oyf_avg / std::pow(numTpls, 1);

	// double oyf_avg = 1.1;
	double xPerc = (totalxLength / numVec);
	double yPerc = (totalyLength / numVec);
	double zPerc = (totalzLength / numVec);

	// cout << " xP: " << xPerc << " yP: " << yPerc << " zP: " << zPerc << " oyf: " << oyf_avg << " ";
	cout << " oyf: " << oyf_avg << " ";

	/*
	cout <<" xL: " << totalxLength << " xP: " << xPerc
		<< " yL: " << totalyLength << " yP: " << yPerc
		<< " zL: " << totalzLength << " zP: " << zPerc
		<< " oyf: " << oyf_avg << endl;

	double curvature = CalculateGaussCurvature(skeleton);
	cout << "Curv: " << curvature << endl;
	*/

	// -1 = undefined,	1 = streamwise,	2 = spanwise,	3 = quasi-streamwise
	//	4 = quasi-spanwise,	5 = vertical,	6 = elevated,	7 = Horseshoe/Hairpin

	cout << " numCells: " << numCells << " ";

	if (xPerc > 0.90)
	{
		cout << "streamwise" << endl;
		std::string type = "streamwise";
		incrementType(vortexCounts, type);
		return 1;
	}
	else if (yPerc > 0.90)
	{
		cout << "spanwise" << endl;
		std::string type = "spanwise";
		incrementType(vortexCounts, type);
		return 2;
	}
	else if (zPerc > 0.90)
	{
		cout << "vertical" << endl;
		std::string type = "vertical";
		incrementType(vortexCounts, type);
		return 5;
	}

	if (oyf_avg < 4)
	{
		if (xPerc < 0.9 && xPerc > yPerc && xPerc > zPerc)
		{
			cout << "quasi-streamwise" << endl;
			std::string type = "quasi-streamwise";
			incrementType(vortexCounts, type);
			return 3;
		}
		else if (yPerc < 0.9 && yPerc > xPerc && yPerc > zPerc)
		{
			cout << "quasi-spanwise" << endl;
			std::string type = "quasi-spanwise";
			incrementType(vortexCounts, type);
			return 4;
		}
		else if (zPerc < 0.9 && zPerc > xPerc && zPerc > yPerc)
		{
			cout << "elevated" << endl;
			std::string type = "elevated";
			incrementType(vortexCounts, type);
			return 6;
		}
	}
	else if (oyf_avg >= 4)
	{
		double curvature = CalculateGaussCurvature(skeleton);
		double gauss_curvature = curvature * (1 - xPerc) * (yPerc) * (zPerc);
		double bbox_length = skeleton->GetLength();

		double 	cor[3], max[3], mid[3], min[3], size[3];
		vtkOBBTree::ComputeOBB(skeleton->GetPoints(), cor, max, mid, min, size);
		/*
		cout << "Corner: " << cor[0] << " " << cor[1] << " " << cor[2] << endl;
		cout << "Max: " << max[0] << " " << max[1] << " " << max[2] << endl;
		cout << "Mid: " << mid[0] << " " << mid[1] << " " << mid[2] << endl;
		cout << "Min: " << min[0] << " " << min[1] << " " << min[2] << endl;
		cout << "Size: " << size[0] << " " << size[1] << " " << size[2] << endl;
		bbox_length = sqrt(pow(sqrt(pow(size[0], 2) + pow(size[1], 2)), 2) + pow(size[2], 2));
		*/
		double p1[3] = { cor[0] + max[0], cor[1] + max[1], cor[2] + max[2] };
		double p2[3] = { cor[0] + mid[0], cor[1] + mid[1], cor[2] + mid[2] };
		double p3[3] = { cor[0] + min[0], cor[1] + min[1], cor[2] + min[2] };
		double a = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
		double b = sqrt(vtkMath::Distance2BetweenPoints(p1, p3));
		double c = sqrt(vtkMath::Distance2BetweenPoints(p2, p3));
		// Pythagorous Theorem
		bbox_length = sqrt(pow(sqrt(pow(a, 2) + pow(b, 2)), 2) + pow(c, 2));

		/*
		if (totallength / len < 0.04)
		{
			cout << "Undefined" << endl;
			return -1;
		}
		*/

		double skelRatio = bbox_length / totallength;
		double newSkelLen = totallength / skelRatio;
		double newSkelRatio = newSkelLen / totallength;
		double normal_length = totallength / std::pow(numVec, 0.25);
		cout << " Gauss Curvature: " << gauss_curvature << " ";
		// cout << " Skel Ratio: " << newSkelRatio << " ";
		// cout << " Normal length: " << normal_length << " ";
		double newCurv = (gauss_curvature * newSkelRatio * normal_length);
		cout << " Curvature: " << newCurv << " ";
		/*
		cout << "Curvature: " << curvature << " " << gauss_curvature << " " << totallength << " "
			<< bbox_length << " " << skelRatio << " " << newSkelRatio << " " << normal_length << " "
			<< newCurv << endl;
		*/

		if (newCurv >= 0.005)
		{
			cout << "Horseshoe/Hairpin" << endl;
			std::string type = "Horseshoe/Hairpin";
			incrementType(vortexCounts, type);
			return 7;
		}
		else if (xPerc < 0.9 && xPerc > yPerc && xPerc > zPerc)
		{
			cout << "quasi-streamwise" << endl;
			std::string type = "quasi-streamwise";
			incrementType(vortexCounts, type);
			return 3;
		}
		else if (yPerc < 0.9 && yPerc > xPerc && yPerc > zPerc)
		{
			cout << "quasi-spanwise" << endl;
			std::string type = "quasi-spanwise";
			incrementType(vortexCounts, type);
			return 7;
		}
		else if (zPerc < 0.9 && zPerc > xPerc && zPerc > yPerc)
		{
			cout << "elevated" << endl;
			std::string type = "elevated";
			incrementType(vortexCounts, type);
			return 6;
		}
		else
		{
			cout << "Undefined" << endl;
			std::string type = "Undefined";
			incrementType(vortexCounts, type);
			return -1;
		}
	}
	else
	{
		cout << "Undefined" << endl;
		std::string type = "Undefined";
		incrementType(vortexCounts, type);
		return -1;
	}
}

int processSkeleton(vtkDataSet* skeleton, const double & zMin, const double & zMax,
	const double & len, std::vector<VortexInfo>& vortexCounts)
{
	vtkNew<vtkNamedColors> colors;

	vtkNew<vtkGeometryFilter> geometry;
	geometry->SetInputData(skeleton);
	geometry->MergingOn();
	geometry->Update();

	vtkNew<vtkCleanPolyData> clean;
	clean->SetInputData(geometry->GetOutput());
	clean->Update();

	vtkNew<vtkStripper> stripper;
	stripper->SetInputData(clean->GetOutput());
	stripper->JoinContiguousSegmentsOn();
	stripper->Update();

	clean->SetInputData(stripper->GetOutput());
	clean->Update();
	vtkPolyData* newSkeleton = clean->GetOutput();
	/*
	int numCells = clean->GetOutput()->GetNumberOfCells();
	for (int i = 0; i < numCells; i++)
	{
		cout << *clean->GetOutput()->GetCell(i) << endl;
	}

	// fixSkeleton(newSkeleton, newSkeleton);
	// cout << *newSkeleton << endl;

	int numCells = skeleton->GetNumberOfCells();
	for (int i = 0; i < numCells; i++)
	{
		cout << *skeleton->GetCell(i) << endl;
	}
	*/

	int skel_category = inspectSkeletonOrientation(newSkeleton, zMin, zMax, len, vortexCounts);
	/*
	if (skel_category != 0)
	{
		return skel_category;
	}

	// Check curvature
	double gauss_curvature = CalculateGaussCurvature(skeleton);
	// cout << "Curvature: " << gauss_curvature << " " << skel_category << endl;
	if (gauss_curvature > 10)
	{
		return 0;
	}
	else
	{
		return -1;
	}


	vtkNew<vtkRenderer> ren1;
	ren1->SetBackground(colors->GetColor3d("White").GetData());

	vtkNew<vtkDataSetMapper> skelMapper;
	skelMapper->SetInputData(skeleton);
	skelMapper->ScalarVisibilityOff();
	skelMapper->Update();
	vtkNew<vtkActor> skelActor;
	skelActor->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());
	skelActor->GetProperty()->SetLineWidth(2.0);
	skelActor->GetProperty()->SetOpacity(1);
	skelActor->SetMapper(skelMapper);
	ren1->AddActor(skelActor);

	ren1->GetActiveCamera()->SetPosition(0, -1, 0);
	ren1->GetActiveCamera()->SetViewUp(0, 0, 1);
	ren1->ResetCamera();

	vtkNew<vtkRenderWindow> ren_win;
	ren_win->SetSize(1200, 600);
	ren_win->SetWindowName("Vortex Lines");
	ren_win->AddRenderer(ren1);

	vtkNew<vtkRenderWindowInteractor> iren;
	vtkNew<vtkInteractorStyleTrackballCamera> irenStyle;
	iren->SetInteractorStyle(irenStyle);
	iren->SetRenderWindow(ren_win);

	ren_win->Render();
	iren->Start();
	*/
	return skel_category;
}


int main(int argc, char* argv[])
{
	vtkSMPTools::SetBackend("STDThread");
	// vtkSMPTools::Initialize(8);
	int num_threads = vtkSMPTools::GetEstimatedNumberOfThreads();
	cout << "Total threads: " << num_threads << endl;

	vtkNew<vtkNamedColors> colors;

	if (argc < 5)
	{
		std::cout << "Please specify the following arguments." << endl;
		std::cout << "1: input skeleton filename." << endl;
		std::cout << "2: input surface filename." << endl;
		std::cout << "3: output vortex type." << endl;
		std::cout << "4: output skeleton filename." << endl;
		std::cout << "5: output surface filename." << endl;
		return EXIT_FAILURE;
	}

	/*
	// Specify the type of vortex statistics to get
	// 0 = all,	1 = streamwise,	2 = spanwise,	3 = quasi-streamwise
	// 4 = quasi-spanwise,	5 = vertical,	6 = elevated,	7 = Omega/Hairpin
	int specifyVortexType = 0;
	// Specify the streamwise direction
	// 0 = X,	1 = Y,	2 = Z
	int specifyStwDir = 0;
	*/

	// std::string datafile = "./skeleton_extraction/fortSkeleton.vtk";
	std::string datafile = argv[1];
	vtkNew<vtkDataSetReader> dataReader;
	dataReader->SetFileName(datafile.c_str());
	dataReader->Update();
	vtkSmartPointer<vtkDataSet> skel_dataset = vtkPolyData::SafeDownCast(dataReader->GetOutput());

	// datafile = "./skeleton_extraction/fortSurface.vtk";
	datafile = argv[2];
	vtkNew<vtkDataSetReader> dataReader2;
	dataReader2->SetFileName(datafile.c_str());
	dataReader2->Update();
	vtkSmartPointer<vtkDataSet> surf_dataset = vtkPolyData::SafeDownCast(dataReader2->GetOutput());
	double dataset_leng = surf_dataset->GetLength();

	double bounds[6];
	skel_dataset->GetBounds(bounds);

	vtkSmartPointer<vtkDataArray> skel_ids_array = skel_dataset->GetCellData()->GetArray("skeletonRegionIds");
	vtkNew<vtkIdList> unique_skel_ids;
	vtkIdType numCells = skel_ids_array->GetNumberOfTuples();

	for (vtkIdType i = 0; i < numCells; i++)
	{
		int id = skel_ids_array->GetTuple1(i);
		unique_skel_ids->InsertUniqueId(id);
	}
	unique_skel_ids->Sort();

	// ############## For Visualization ####################
	vtkNew<vtkRenderer> ren1;
	vtkNew<vtkRenderWindow> ren_win;
	ren_win->SetSize(1200, 1000);
	ren_win->SetWindowName("Vortex Lines");
	ren_win->AddRenderer(ren1);

	vtkNew<vtkDataSetMapper> dataMapper;
	dataMapper->SetInputData(skel_dataset);
	dataMapper->ScalarVisibilityOff();
	dataMapper->Update();
	vtkNew<vtkActor> dataActor;
	dataActor->GetProperty()->SetColor(colors->GetColor3d("Red").GetData());
	dataActor->GetProperty()->SetLineWidth(2.0);
	dataActor->GetProperty()->SetOpacity(1);
	dataActor->SetMapper(dataMapper);
	ren1->AddActor(dataActor);

	ren1->GetActiveCamera()->SetPosition(0, -1, 0);
	ren1->GetActiveCamera()->SetViewUp(0, 0, 1);
	ren1->ResetCamera();
	ren1->GetActiveCamera()->Azimuth(-30);
	ren1->GetActiveCamera()->Elevation(15);
	ren1->GetActiveCamera()->Zoom(1.65);

	vtkNew<vtkRenderWindowInteractor> iren;
	vtkNew<vtkInteractorStyleTrackballCamera> irenStyle;
	iren->SetInteractorStyle(irenStyle);
	iren->SetRenderWindow(ren_win);

	// ############## End For Visualization ####################

	std::vector<VortexInfo> vortexCounts;

	int skel_category;
	std::string arr_name = "RegionIds";
	std::vector<Region<vtkUnstructuredGrid>> skelDatasets;
	std::vector<Region<vtkUnstructuredGrid>> surfDatasets;
	// First divide the ids into blocks and get the datasets into the vector
	if (unique_skel_ids->GetNumberOfIds() > 5 * num_threads)
	{
		int stepSize = 5 * num_threads; // numRegions / numThreads;
		int startRegionIdx = 0, endRegionIdx = 0;
		int startRegionId = 0, endRegionId = 0;

		vtkNew<vtkThreshold> threshold;
		threshold->SetInputData(skel_dataset);
		threshold->SetInputArrayToProcess(0, 0, 0, 1, "skeletonRegionIds");

		vtkNew<vtkThreshold> threshold2;
		threshold2->SetInputData(surf_dataset);
		threshold2->SetInputArrayToProcess(0, 0, 0, 1, arr_name.c_str());

		int stepId = -1;
		bool toBreak = false;
		while (true)
		{
			stepId++;
			startRegionIdx = endRegionIdx;
			endRegionIdx = startRegionIdx + stepSize;
			startRegionId = unique_skel_ids->GetId(startRegionIdx);
			endRegionId = unique_skel_ids->GetId(endRegionIdx - 1);

			if (endRegionIdx >= unique_skel_ids->GetNumberOfIds())
			{
				endRegionIdx = unique_skel_ids->GetNumberOfIds();
				endRegionId = unique_skel_ids->GetId(endRegionIdx - 1);
				toBreak = true;
			}

			cout << stepId << " " << startRegionId << " " << endRegionId << endl;

			threshold->SetLowerThreshold(startRegionId);
			threshold->SetUpperThreshold(endRegionId);
			threshold->Update();

			threshold2->SetLowerThreshold(startRegionId);
			threshold2->SetUpperThreshold(endRegionId);
			threshold2->Update();

			// vtkSmartPointer<vtkUnstructuredGrid> tempGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
			// tempGrid->DeepCopy(thresholdFil->GetOutput());
			// cout << "Creating new region object..." << endl;

			// Store this region into a vector for later retrival
			Region<vtkUnstructuredGrid> region1(stepId, threshold->GetOutput());
			skelDatasets.push_back(region1);

			Region<vtkUnstructuredGrid> region2(stepId, threshold2->GetOutput());
			surfDatasets.push_back(region2);
			// cout << "Adding new region object to the list..." << endl;
			if (toBreak)
			{
				break;
			}
		}
	}
	else
	{
		vtkSmartPointer<vtkUnstructuredGrid> tempGrid1 = vtkSmartPointer<vtkUnstructuredGrid>::New();
		tempGrid1->DeepCopy(skel_dataset);
		Region<vtkUnstructuredGrid> region1(0, tempGrid1);
		skelDatasets.push_back(region1);

		vtkSmartPointer<vtkUnstructuredGrid> tempGrid2 = vtkSmartPointer<vtkUnstructuredGrid>::New();
		tempGrid2->DeepCopy(surf_dataset);
		Region<vtkUnstructuredGrid> region2(0, tempGrid2);
		surfDatasets.push_back(region2);
	}

	vtkSmartPointer<vtkUnstructuredGrid> full_skeleton_data = nullptr;
	vtkSmartPointer<vtkUnstructuredGrid> full_surface_data = nullptr;
	for (int datasetId = 0; datasetId < skelDatasets.size(); datasetId++)
	{
		vtkSmartPointer<vtkUnstructuredGrid> outDataset = nullptr;
		vtkSmartPointer<vtkUnstructuredGrid> outSurfDataset = nullptr;

		vtkUnstructuredGrid* tempGrid0 = skelDatasets.at(datasetId).RegionData;
		vtkSmartPointer<vtkThreshold> threshold = vtkSmartPointer<vtkThreshold>::New();
		threshold->SetInputData(tempGrid0);
		threshold->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, "skeletonRegionIds");

		vtkUnstructuredGrid* tempGrid1 = surfDatasets.at(datasetId).RegionData;
		vtkSmartPointer<vtkThreshold> threshold2 = vtkSmartPointer<vtkThreshold>::New();
		threshold2->SetInputData(tempGrid1);
		threshold2->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, arr_name.c_str());

		double idsRange[2];
		tempGrid0->GetCellData()->GetArray("skeletonRegionIds")->GetRange(idsRange);

		for (int id = idsRange[0]; id <= idsRange[1]; id++)
		{
			vtkIdType idloc = unique_skel_ids->FindIdLocation(id);
			if (idloc == -1)
			{
				continue;
			}
			cout << "Skel Num: " << id << " ";
			threshold->SetLowerThreshold(id);
			threshold->SetUpperThreshold(id);
			threshold->Update();
			vtkSmartPointer<vtkDataSet> skeleton_data = threshold->GetOutput();

			if (skeleton_data->GetNumberOfPoints() == 0)
			{
				continue;
			}
			skel_category = processSkeleton(skeleton_data, bounds[4], bounds[5],
				dataset_leng, vortexCounts);
			if (skel_category == -1)
			{
				continue;
			}

			threshold2->SetLowerThreshold(id);
			threshold2->SetUpperThreshold(id);
			threshold2->Update();
			vtkSmartPointer<vtkDataSet> surface_data = threshold2->GetOutput();

			/**/
			if (std::stoi(argv[3]) == 0)
			{
				// Do nothing
			}
			else if (skel_category != std::stoi(argv[3]))
			{
				continue;
			}

			if (outDataset == nullptr)
			{
				outDataset = vtkSmartPointer<vtkUnstructuredGrid>::New();
				outDataset->DeepCopy(skeleton_data);

				outSurfDataset = vtkSmartPointer<vtkUnstructuredGrid>::New();
				outSurfDataset->DeepCopy(surface_data);
			}
			else
			{
				vtkSmartPointer<vtkAppendFilter> appendFilter = vtkSmartPointer<vtkAppendFilter>::New();
				appendFilter->AddInputData(outDataset);
				appendFilter->AddInputData(skeleton_data);
				appendFilter->Update();
				outDataset->DeepCopy(appendFilter->GetOutput());

				appendFilter = vtkSmartPointer<vtkAppendFilter>::New();
				appendFilter->AddInputData(outSurfDataset);
				appendFilter->AddInputData(surface_data);
				appendFilter->Update();
				outSurfDataset->DeepCopy(appendFilter->GetOutput());
			}

			/*
			// ############## Visualization Start ####################
			// All the stuff for the first renderer.

			ren1->SetBackground(colors->GetColor3d("White").GetData());
			ren1->SetViewport(0, 0, 1, 0.5);
			// vtkCamera* camera = ren1->GetActiveCamera();

			vtkNew<vtkRenderer> ren2;
			ren2->SetBackground(colors->GetColor3d("White").GetData());
			ren2->SetViewport(0, 0.5, 1, 1);
			// ren2->SetActiveCamera(camera);

			vtkNew<vtkDataSetMapper> vortexMapper;
			vortexMapper->SetInputData(skeleton_data);
			vortexMapper->ScalarVisibilityOff();
			vortexMapper->Update();
			vtkNew<vtkActor> vortexActor;
			vortexActor->GetProperty()->SetColor(colors->GetColor3d("Blue").GetData());
			vortexActor->GetProperty()->SetLineWidth(2.0);
			vortexActor->GetProperty()->SetOpacity(1);
			vortexActor->SetMapper(vortexMapper);

			vtkNew<vtkDataSetMapper> vortexSurfMapper;
			vortexSurfMapper->SetInputData(surface_data);
			vortexSurfMapper->ScalarVisibilityOff();
			vortexSurfMapper->Update();
			vtkNew<vtkActor> vortexSurfActor;
			vortexSurfActor->GetProperty()->SetColor(colors->GetColor3d("SkyBlue").GetData());
			vortexSurfActor->GetProperty()->SetLineWidth(2.0);
			vortexSurfActor->GetProperty()->SetOpacity(1);
			vortexSurfActor->SetMapper(vortexSurfMapper);

			vtkNew<vtkOutlineFilter> outline;
			outline->SetInputData(surface_data);
			outline->Update();
			vtkNew<vtkDataSetMapper> bboxMapper;
			bboxMapper->SetInputData(outline->GetOutput());
			bboxMapper->ScalarVisibilityOff();
			bboxMapper->Update();
			vtkNew<vtkActor> bboxActor;
			bboxActor->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());
			bboxActor->GetProperty()->SetLineWidth(2.0);
			bboxActor->GetProperty()->SetOpacity(0.5);
			bboxActor->SetMapper(bboxMapper);

			ren1->AddActor(vortexActor);
			ren1->AddActor(vortexSurfActor);
			ren1->AddActor(bboxActor);

			// All the stuff for the 2nd renderer

			vtkNew<vtkDataSetMapper> dataMapper2;
			dataMapper2->SetInputData(skeleton_data);
			dataMapper2->ScalarVisibilityOff();
			dataMapper2->Update();
			vtkNew<vtkActor> dataActor2;
			dataActor2->GetProperty()->SetColor(colors->GetColor3d("Black").GetData());
			dataActor2->GetProperty()->SetLineWidth(2.0);
			dataActor2->GetProperty()->SetOpacity(1);
			dataActor2->SetMapper(dataMapper2);

			vtkNew<vtkDataSetMapper> surfDataMapper;
			surfDataMapper->SetInputData(surface_data);
			surfDataMapper->ScalarVisibilityOff();
			surfDataMapper->Update();
			vtkNew<vtkActor> surfDataActor;
			surfDataActor->GetProperty()->SetColor(colors->GetColor3d("SkyBlue").GetData());
			surfDataActor->GetProperty()->SetOpacity(0.5);
			surfDataActor->SetMapper(surfDataMapper);

			vtkColor3d axesActorColor = colors->GetColor3d("Black");
			vtkNew<vtkCubeAxesActor> cubeAxesActor;
			cubeAxesActor->SetUseTextActor3D(0);
			cubeAxesActor->SetBounds(outline->GetOutput()->GetBounds());
			cubeAxesActor->SetLabelOffset(5);
			cubeAxesActor->SetTitleOffset(5);
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
			cubeAxesActor->SetCamera(ren2->GetActiveCamera());

			ren2->AddActor(dataActor2);
			ren2->AddActor(surfDataActor);
			ren2->AddActor(cubeAxesActor);

			ren2->GetActiveCamera()->SetPosition(0, -1, 0);
			ren2->GetActiveCamera()->SetViewUp(0, 0, 1);
			ren2->ResetCamera();
			ren2->GetActiveCamera()->Azimuth(-30);
			ren2->GetActiveCamera()->Elevation(15);
			ren2->GetActiveCamera()->Zoom(1);

			ren_win->AddRenderer(ren2);
			ren_win->Render();
			iren->Start();
			ren_win->RemoveRenderer(ren2);

			ren1->RemoveActor(vortexActor);
			ren1->RemoveActor(vortexSurfActor);
			ren1->RemoveActor(bboxActor);

			ren2->RemoveActor(dataActor2);
			ren2->RemoveActor(surfDataActor);
			ren2->RemoveActor(cubeAxesActor);
			*/
			// ############## Visualization End ####################
		}

		if (outDataset == nullptr || outSurfDataset == nullptr)
		{
			continue;
		}

		if (full_skeleton_data == nullptr)
		{
			full_skeleton_data = vtkSmartPointer<vtkUnstructuredGrid>::New();
			full_skeleton_data->DeepCopy(outDataset);
			full_surface_data = vtkSmartPointer<vtkUnstructuredGrid>::New();
			full_surface_data->DeepCopy(outSurfDataset);
		}
		else
		{
			vtkSmartPointer<vtkAppendFilter> appendFilter = vtkSmartPointer<vtkAppendFilter>::New();
			appendFilter->AddInputData(full_skeleton_data);
			appendFilter->AddInputData(outDataset);
			appendFilter->Update();
			full_skeleton_data->DeepCopy(appendFilter->GetOutput());

			appendFilter = vtkSmartPointer<vtkAppendFilter>::New();
			appendFilter->AddInputData(full_surface_data);
			appendFilter->AddInputData(outSurfDataset);
			appendFilter->Update();
			full_surface_data->DeepCopy(appendFilter->GetOutput());
		}
	}

	cout << endl << "#################################" << endl;
	cout << "#################################" << endl;
	cout << "Total: " << unique_skel_ids->GetNumberOfIds() << endl;
	for (VortexInfo vortexCount : vortexCounts)
	{
		cout << vortexCount.getType() << ": " << vortexCount.getCount() << endl;
	}

	vtkNew<vtkDataSetWriter> writer;
	writer->SetInputData(full_skeleton_data);
	// std::string temp = "./vortex_statistics/fortSkel.vtk";
	std::string temp = argv[4];
	writer->SetFileName(temp.c_str());
	writer->SetFileTypeToBinary();
	writer->Write();

	vtkNew<vtkDataSetWriter> writer2;
	writer2->SetInputData(full_surface_data);
	// temp = "./vortex_statistics/fortSurf.vtk";
	temp = argv[5];
	writer2->SetFileName(temp.c_str());
	writer2->SetFileTypeToBinary();
	writer2->Write();

	return EXIT_SUCCESS;
}