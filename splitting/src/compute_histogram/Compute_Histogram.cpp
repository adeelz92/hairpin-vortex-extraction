// In this version, we extract lambda2 steps for tree construction using histogram strategy.
// Compute_Histogram.exe "Path\To\FullDataVTKfile.vtk" "Path\To\Dataset_steps.csv" 89

#include <vtkAxes.h>
#include <vtkDataSetReader.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

struct Bin {
	int count;
	float lowerRange;
	float upperRange;
};

struct Row {
	int index;
	float lowerRange;
	float upperRange;
	float percentage;
};

std::vector<Bin> computeHistogram(std::vector<float>& values, int& numValues, int& numBins,
	float& minValue, float& maxValue) {
	// Create a vector to store each bin's count and range
	std::vector<Bin> histogram(numBins);

	// Compute the range and bin width
	float range = maxValue - minValue;
	float binWidth = range / numBins;

	// Initialize each bin's count and range
	for (int i = 0; i < numBins; i++) {
		histogram[i].count = 0;
		histogram[i].lowerRange = minValue + (i * binWidth);
		histogram[i].upperRange = std::min(minValue + ((i + 1) * binWidth), maxValue);
	}

	float v;
	// Compute the bin index for each value and increment the corresponding bin count
	for (int i = 0; i < numValues; i++) {
		v = values[i];
		if (v < minValue || v > maxValue)
		{
			continue;
		}
		int binIndex = static_cast<int>((values[i] - minValue) / binWidth);

		// Increment the bin count, ensuring it stays within the valid range
		histogram[std::min(binIndex, numBins - 1)].count++;
	}

	return histogram;
}

void writeRowsToFile(const std::vector<Row>& rows, const std::string& filePath) {
	std::ofstream file(filePath);
	if (!file.is_open()) {
		std::cout << "Failed to open file: " << filePath << std::endl;
		return;
	}

	for (const auto& row : rows) {
		file << row.index << ", " << row.lowerRange << ", " << row.upperRange << ", " << row.percentage << std::endl;
		std::cout << row.index << ", " << row.lowerRange << ", " << row.upperRange << ", " << row.percentage << std::endl;
	}

	file.close();
}

void mainFunction(const std::string& datasetPath, const std::string& stepsSavePath, const int& percentile)
{
	vtkSmartPointer<vtkDataSetReader> reader = vtkSmartPointer<vtkDataSetReader>::New();
	reader->SetFileName(datasetPath.c_str());
	reader->Update();

	vtkSmartPointer<vtkDataSet> dataset = vtkDataSet::SafeDownCast(reader->GetOutput());

	vtkSmartPointer<vtkDataArray> lambda2Arr = dataset->GetPointData()->GetArray("lambda2-vtk");

	std::vector<float> negLambda2;
	vtkIdType totalValues = lambda2Arr->GetNumberOfTuples();
	for (vtkIdType i = 0; i < totalValues; i++)
	{
		float value = lambda2Arr->GetTuple1(i);
		if (value < 0)
		{
			negLambda2.push_back(value);
		}
	}

	int numValues = negLambda2.size();
	int numBins = 100;
	int prevCount = 0;
	// Find the minimum and maximum values in the array
	float minValue = negLambda2[0];
	float maxValue = negLambda2[0];
	for (int i = 1; i < numValues; i++) {
		if (negLambda2[i] < minValue) {
			minValue = negLambda2[i];
		}
		if (negLambda2[i] > maxValue) {
			maxValue = negLambda2[i];
		}
	}
	float minValueBack = minValue;
	std::vector<Bin> histogram;
	while (true) {

		std::vector<float> newBins;
		std::vector<int> newCounts;

		// Compute histogram
		histogram = computeHistogram(negLambda2, numValues, numBins, minValue, maxValue);

		// Find bins with at least 0.01% of values
		for (int i = 0; i < numBins; i++) {
			float perc = (static_cast<float>(histogram[i].count) / numValues) * 100;
			if (perc > 0.1) {
				newBins.push_back(histogram[i].lowerRange);
				newCounts.push_back(histogram[i].count);
			}
		}
		newBins.push_back(histogram[numBins - 1].upperRange);


		// Update variables for next iteration
		minValue = newBins[0];
		maxValue = newBins[newBins.size() - 1];

		// Early stopping conditions
		float lastPerc = (static_cast<float>(histogram[numBins - 1].count) / numValues) * 100;
		float secLastPerc = (static_cast<float>(histogram[numBins - 2].count) / numValues) * 100;

		if ((lastPerc < 30 && lastPerc - secLastPerc < 20) || prevCount == histogram[numBins - 1].count) {
			break;
		}
		prevCount = histogram[numBins - 1].count;
	}

	/*
	// Output the bin count and range
	for (int i = 0; i < numBins; i++) {
		std::cout << "Bin " << i << ": Count = " << histogram[i].count
			<< ", Range = [" << histogram[i].lowerRange << ", " << histogram[i].upperRange << "]" << std::endl;
	}
	*/
	minValue = minValueBack;
	maxValue = histogram[percentile].lowerRange;
	cout << minValue << " " << maxValue << endl;

	while (true) {
		histogram = computeHistogram(negLambda2, numValues, numBins, minValue, maxValue);
		float perc = (static_cast<float>(histogram[numBins - 1].count) / numValues) * 100;
		if (perc < 10) {
			break;
		}
		else {
			numBins *= 2;
		}
	}

	std::vector<float> finalBins;
	std::vector<int> finalCounts;

	// Iterate over histogram bins
	for (int i = 0; i < histogram.size(); i++) {
		float perc = (static_cast<float>(histogram[i].count) / numValues) * 100;

		if (perc <= 0) {
			continue;
		}

		finalBins.insert(finalBins.begin(), histogram[i].lowerRange);
		finalCounts.insert(finalCounts.begin(), histogram[i].count);
	}

	finalBins.insert(finalBins.begin(), histogram[histogram.size() - 1].upperRange);

	std::cout << "Length: " << finalBins.size() << std::endl;

	std::vector<Row> rows;

	// Insert the first element
	Row firstRow = { 0, finalBins[0], finalBins[1], (finalCounts[0] / static_cast<float>(negLambda2.size())) * 100 };
	rows.insert(rows.begin(), firstRow);

	int prevSum = 1;
	int currSum = 1;
	while (currSum < finalCounts.size()) {
		float perc = (finalCounts[currSum] / static_cast<float>(negLambda2.size())) * 100;
		Row row = { currSum, finalBins[currSum], finalBins[currSum + 1], perc };
		rows.insert(rows.begin(), row);

		int temp = currSum;
		currSum = prevSum + currSum;
		prevSum = temp;
	}

	writeRowsToFile(rows, stepsSavePath);

}

int main(int argc, char* argv[])
{
	if (argc < 4)
	{
		std::cerr << "Please specify input and output files." << endl;
		exit(1);
	}
	std::string filepath = argv[1];
	std::string savePath2 = argv[2];
	const int percentile = std::stoi(argv[3]);
	// std::string filepath = "D:\\Adeel\\studies\\UH\\Data_Visualization\\Assignments(C++)\\tutorials\\vortex_verification\\build\\Release\\all_datasets_full\\fortFull_bin.vtk";
	// std::string savePath2 = "D:\\Adeel\\studies\\UH\\Data_Visualization\\Assignments(C++)\\tutorials\\vortex_verification\\build\\Release\\initial_value_problem\\fortFull_bin.csv";
	std::vector<std::string> filenames;

	mainFunction(filepath, savePath2, percentile);

	return 0;
}
