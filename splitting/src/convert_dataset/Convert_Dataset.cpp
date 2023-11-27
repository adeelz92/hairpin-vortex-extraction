#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <vtkStructuredGrid.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkSMPTools.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkDataSetWriter.h>
#include <chrono>
#include <cctype>

// Convert_Dataset.exe "Path\To\RawFile" "Path\To\VTKfile.vtk" xDim yDim zDim
int main(int argc, char* argv[])
{
	vtkSMPTools::SetBackend("STDThread");
	// vtkSMPTools::Initialize(2);
	int num_threads = vtkSMPTools::GetEstimatedNumberOfThreads();
	cout << "Total threads: " << num_threads << endl;

	// cout << argv[1] << endl << argv[2] << endl;
	// std::string src_file_path = "D:\\Adeel\\studies\\UH\\Research\\project-data\\fort.80110010";
	std::ifstream src_file(argv[1]);
	std::vector<std::vector<std::string>> lines;
	std::string line;

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

	cout << "Reading data..." << endl;
	bool containsAlphabets;
	while (std::getline(src_file, line)) {
		std::istringstream iss(line);
		std::vector<std::string> components;
		std::string component;
		containsAlphabets = false;
		while (iss >> component) {

			for (char c : component) {
				if (std::isalpha(c) && c != 'e') // Check for scientific notation
				{
					containsAlphabets = true;
					break; // No need to continue checking
				}
			}

			if (containsAlphabets == true)
			{
				break;
			}

			components.push_back(component);
		}

		if (containsAlphabets == true) {
			continue;
		}

		lines.push_back(components);
	}

	// int num_x = 384;
	int num_x = std::stoi(argv[3]);
	// int num_y = 384;
	int num_y = std::stoi(argv[4]);
	// int num_z = 193;
	int num_z = std::stoi(argv[5]);
	int numTuples = num_x * num_y * num_z;

	cout << "Total valid lines: " << lines.size() << " Total tuples: " << numTuples << endl;

	if (lines.size() != numTuples)
	{
		std::cerr << "Total data lines in the file should be equal to total number of points." << endl;
		exit(1);
	}

	cout << "Converting data..." << endl;
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(numTuples);
	vtkSmartPointer<vtkFloatArray> velocity_arr = vtkSmartPointer<vtkFloatArray>::New();
	velocity_arr->SetName("velocity");
	velocity_arr->SetNumberOfComponents(3);
	velocity_arr->SetNumberOfTuples(numTuples);
	vtkSmartPointer<vtkFloatArray> velocity_uf_arr = vtkSmartPointer<vtkFloatArray>::New();
	velocity_uf_arr->SetName("velocity-uf");
	velocity_uf_arr->SetNumberOfComponents(3);
	velocity_uf_arr->SetNumberOfTuples(numTuples);
	vtkSmartPointer<vtkFloatArray> vorticity_arr = vtkSmartPointer<vtkFloatArray>::New();
	vorticity_arr->SetName("vorticity");
	vorticity_arr->SetNumberOfComponents(3);
	vorticity_arr->SetNumberOfTuples(numTuples);
	vtkSmartPointer<vtkFloatArray> vorticity_oyf_arr = vtkSmartPointer<vtkFloatArray>::New();
	vorticity_oyf_arr->SetName("vorticity-oyf");
	vorticity_oyf_arr->SetNumberOfComponents(3);
	vorticity_oyf_arr->SetNumberOfTuples(numTuples);
	vtkSmartPointer<vtkFloatArray> oyf_arr = vtkSmartPointer<vtkFloatArray>::New();
	oyf_arr->SetName("oyf");
	oyf_arr->SetNumberOfComponents(1);
	oyf_arr->SetNumberOfTuples(numTuples);
	vtkSmartPointer<vtkFloatArray> lambda2_arr = vtkSmartPointer<vtkFloatArray>::New();
	lambda2_arr->SetName("lambda2");
	lambda2_arr->SetNumberOfComponents(1);
	lambda2_arr->SetNumberOfTuples(numTuples);
	cout << "Allocate finished... " << endl;
	vtkSMPTools::For(0, lines.size(), [&](vtkIdType tupleId, vtkIdType endId) {
		auto& Velocity_arr = vtk::DataArrayTupleRange<3>(velocity_arr);
		auto& Velocity_u_arr = vtk::DataArrayTupleRange<3>(velocity_uf_arr);
		auto& Vorticity_arr = vtk::DataArrayTupleRange<3>(vorticity_arr);
		auto& Vorticity_oyf_arr = vtk::DataArrayTupleRange<3>(vorticity_oyf_arr);
		auto& Oyf_arr = vtk::DataArrayTupleRange<1>(oyf_arr);
		auto& Lambda2_arr = vtk::DataArrayTupleRange<1>(lambda2_arr);
		for (; tupleId < endId; ++tupleId) {
			const auto& components = lines[tupleId];
			float x = std::stof(components[0]);
			float y = std::stof(components[1]);
			float z = std::stof(components[2]);
			float u = std::stof(components[3]);
			float uf = std::stof(components[4]);
			float v = std::stof(components[5]);
			float w = std::stof(components[6]);
			float lambda2 = std::stof(components[7]);
			float ox = std::stof(components[8]);
			float oy = std::stof(components[9]);
			float oyf = std::stof(components[10]);
			float oz = std::stof(components[11]);

			float point[3] = { x, y, z };
			points->InsertPoint(tupleId, point);

			Velocity_arr[tupleId][0] = u;
			Velocity_arr[tupleId][1] = v;
			Velocity_arr[tupleId][2] = w;

			Velocity_u_arr[tupleId][0] = uf;
			Velocity_u_arr[tupleId][1] = v;
			Velocity_u_arr[tupleId][2] = w;

			Vorticity_arr[tupleId][0] = ox;
			Vorticity_arr[tupleId][1] = oy;
			Vorticity_arr[tupleId][2] = oz;

			Vorticity_oyf_arr[tupleId][0] = ox;
			Vorticity_oyf_arr[tupleId][1] = oyf;
			Vorticity_oyf_arr[tupleId][2] = oz;

			Oyf_arr[tupleId][0] = oyf;

			Lambda2_arr[tupleId][0] = lambda2;
		}
	});
	cout << "Insertion finished... " << endl;

	vtkSmartPointer<vtkStructuredGrid> dataset = vtkSmartPointer<vtkStructuredGrid>::New();
	dataset->SetDimensions(num_x, num_y, num_z);
	dataset->SetPoints(points);
	dataset->GetPointData()->AddArray(velocity_arr);
	dataset->GetPointData()->AddArray(velocity_uf_arr);
	dataset->GetPointData()->AddArray(vorticity_arr);
	dataset->GetPointData()->AddArray(vorticity_oyf_arr);
	dataset->GetPointData()->AddArray(oyf_arr);
	dataset->GetPointData()->AddArray(lambda2_arr);
	dataset->Squeeze();

	cout << "Writing data to the vtk file..." << endl;
	// std::string dst_file_path = "D:\\Adeel\\studies\\UH\\Research\\project-data\\fort180_bin.vtk";
	vtkSmartPointer<vtkDataSetWriter> writer = vtkSmartPointer<vtkDataSetWriter>::New();
	writer->SetFileName(argv[2]);
	writer->SetInputData(dataset);
	writer->SetFileTypeToBinary();
	writer->Write();

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

	return 0;
}