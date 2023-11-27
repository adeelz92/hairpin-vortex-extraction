// In this version of checkpoint, we use the velocity to compute jacobian,
// acceleration, q, delta, lambda2, lambdaci, vorticity and oyf values
// Create_Full_Dataset.exe "Path\To\VTKfile.vtk" "Path\To\FullDataVTKfile.vtk"
#include <vtkArrayDispatch.h>
#include <vtkSMPTools.h>
#include <vtkDataSetReader.h>
#include <vtkGradientFilter.h>
#include <vtkDataObject.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSet.h>
#include <vtkStructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkGenerateGlobalIds.h>
#include <vtkDataSetWriter.h>
#include <vtk_eigen.h>
#include VTK_EIGEN(Eigenvalues)
#include VTK_EIGEN(Geometry)

#include <array>
#include <iostream>
#include <conio.h>

namespace
{
	// Computes A*b = x given a 3-matrix A and a 3-vector b.
	template <typename AArrayType, typename BArrayType, typename XArrayType>
	class MatrixVectorMultiplyFunctor
	{
		AArrayType* AArray;
		BArrayType* BArray;
		XArrayType* XArray;

	public:
		MatrixVectorMultiplyFunctor(AArrayType* aArray, BArrayType* bArray, XArrayType* xArray)
			: AArray(aArray)
			, BArray(bArray)
			, XArray(xArray)
		{
		}

		void operator()(vtkIdType begin, vtkIdType end)
		{
			const auto aRange = vtk::DataArrayTupleRange<9>(this->AArray, begin, end);
			const auto bRange = vtk::DataArrayTupleRange<3>(this->BArray, begin, end);
			auto xRange = vtk::DataArrayTupleRange<3>(this->XArray, begin, end);

			auto a = aRange.cbegin();
			auto b = bRange.cbegin();
			auto x = xRange.begin();

			for (; a != aRange.cend(); ++a, ++b, ++x)
			{
				for (vtkIdType i = 0; i < 3; ++i)
				{
					(*x)[i] =
						((*a)[0 + i * 3] * (*b)[0] + (*a)[1 + i * 3] * (*b)[1] + (*a)[2 + i * 3] * (*b)[2]);
				}
			}
		}
	};

	struct MatrixVectorMultiplyWorker
	{
		template <typename AArrayType, typename BArrayType, typename XArrayType>
		void operator()(AArrayType* aArray, BArrayType* bArray, XArrayType* xArray)
		{
			MatrixVectorMultiplyFunctor<AArrayType, BArrayType, XArrayType> functor(aArray, bArray, xArray);
			vtkSMPTools::For(0, xArray->GetNumberOfTuples(), functor);
		}
	};


	void computeVortexCriteria(const float s[9], const float omega[9], float vortexCriteria[4])
	{

		Eigen::Matrix<double, 3, 3> S, Omega, J;
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				const double& s_ij = s[3 * i + j];
				const double& omega_ij = omega[3 * i + j];
				S(i, j) = s_ij;
				Omega(i, j) = omega_ij;
				J(i, j) = (s_ij + omega_ij) / 2.;
			}
		}

		// The Q-criterion is defined as
		// Q = \frac{1}{2} \left[ | \Omega |^2 - | S |^2 \right] > 0
		float& Q = vortexCriteria[0];
		Q = (Omega.operatorNorm() - S.operatorNorm()) / 2.;

		// The delta-criterion is defined as
		// \Delta = \left( \frac{Q}{3} \right)^3 + \left( \frac{\det J}{2} \right)^2 > 0
		float& delta = vortexCriteria[1];
		delta = std::pow(Q / 3., 3) + std::pow(J.determinant() / 2., 2);

		// The lambda_2-criterion is defined as
		// \lambda_2 \left( S^2 + \Omega^2 \right) < 0
		// where $\lambda_2$ is the intermediate eigenvalue
		float& lambda_2 = vortexCriteria[2];
		{
			Eigen::Matrix<double, 3, 3> A = S * S + Omega * Omega;
			auto Eigenvalues = A.eigenvalues();
			// Matrix A is symmetric, so its eigenvalues are all real.
			std::array<double, 3> eigenvalues = { Eigenvalues[0].real(), Eigenvalues[1].real(),
			  Eigenvalues[2].real() };
			std::nth_element(eigenvalues.begin(), eigenvalues.begin() + 1, eigenvalues.end());
			lambda_2 = eigenvalues[1];
		}

		// The lambda_ci-criterion is defined as the imaginary component of the
	  // eigenvalues of the complex conjugate pair of eigenvalues of J
		float& lambda_ci = vortexCriteria[3];
		{
			Eigen::EigenSolver<Eigen::Matrix<double, 3, 3>> eigensolver(J);
			auto eigenvalues = eigensolver.eigenvalues();

			if (std::abs(eigenvalues[0].imag()) > VTK_DBL_EPSILON)
			{
				if ((std::abs(eigenvalues[0].real() - eigenvalues[1].real()) < VTK_DBL_EPSILON &&
					std::abs(eigenvalues[0].imag() + eigenvalues[1].imag()) < VTK_DBL_EPSILON) ||
					(std::abs(eigenvalues[0].real() - eigenvalues[2].real()) < VTK_DBL_EPSILON &&
						std::abs(eigenvalues[0].imag() + eigenvalues[2].imag()) < VTK_DBL_EPSILON))
				{
					lambda_ci = std::abs(eigenvalues[0].imag());
				}
			}
			else if (std::abs(eigenvalues[1].imag()) > VTK_DBL_EPSILON)
			{
				if (std::abs(eigenvalues[1].real() - eigenvalues[2].real()) < VTK_DBL_EPSILON &&
					std::abs(eigenvalues[1].imag() + eigenvalues[2].imag()) < VTK_DBL_EPSILON)

				{
					lambda_ci = std::abs(eigenvalues[1].imag());
				}
			}
		}
	}


	template <typename JacobianArrayType, typename QArrayType, typename DeltaArrayType,
		typename Lambda2ArrayType, typename LambdaciArrayType>
		class ComputeCriteriaFunctor
	{
		JacobianArrayType* JacobianArray;
		QArrayType* QArray;
		DeltaArrayType* DeltaArray;
		Lambda2ArrayType* Lambda2Array;
		LambdaciArrayType* LambdaciArray;

	public:
		ComputeCriteriaFunctor(JacobianArrayType* jacobianArray, QArrayType* qArray,
			DeltaArrayType* deltaArray, Lambda2ArrayType* lambda2Array, LambdaciArrayType* lambdaciArray)
		{
			JacobianArray = jacobianArray;
			QArray = qArray;
			DeltaArray = deltaArray;
			Lambda2Array = lambda2Array;
			LambdaciArray = lambdaciArray;
		}

		void operator() (vtkIdType begin, vtkIdType end)
		{
			const auto jacobianRange = vtk::DataArrayTupleRange<9>(this->JacobianArray, begin, end);
			auto qRange = vtk::DataArrayTupleRange<1>(this->QArray, begin, end);
			auto deltaRange = vtk::DataArrayTupleRange<1>(this->DeltaArray, begin, end);
			auto lambda2Range = vtk::DataArrayTupleRange<1>(this->Lambda2Array, begin, end);
			auto lambdaciRange = vtk::DataArrayTupleRange<1>(this->LambdaciArray, begin, end);

			auto j = jacobianRange.cbegin();
			auto q = qRange.begin();
			auto delta = deltaRange.begin();
			auto lambda2 = lambda2Range.begin();
			auto lambdaci = lambdaciRange.begin();

			for (; j != jacobianRange.cend(); ++j, ++q, ++delta, ++lambda2, ++lambdaci)
			{
				std::array<float, 4> vortexCriteria;
				float S[9];
				float Omega[9];
				static const std::array<typename decltype(j)::value_type::size_type, 9> idxT = { 0, 3, 6, 1,
					4, 7, 2, 5, 8 };
				for (int i = 0; i < 9; i++)
				{
					double j_i = (*j)[i];
					double jt_i = (*j)[idxT[i]];

					S[i] = (j_i + jt_i) / 2.;
					Omega[i] = (j_i - jt_i) / 2.;
				}

				computeVortexCriteria(S, Omega, vortexCriteria.data());
				(*q)[0] = vortexCriteria[0];
				(*delta)[0] = vortexCriteria[1];
				(*lambda2)[0] = vortexCriteria[2];
				(*lambdaci)[0] = vortexCriteria[3];
			}
		}
	};

	struct ComputeCriteriaWorker
	{
		template <typename JacobianArrayType, typename QArrayType, typename DeltaArrayType,
			typename Lambda2ArrayType, typename LambdaciArrayType>
			void operator()(JacobianArrayType* jacobianArray, QArrayType* qArray,
				DeltaArrayType* deltaArray, Lambda2ArrayType* lambda2Array, LambdaciArrayType* lambdaciArray)
		{
			ComputeCriteriaFunctor<JacobianArrayType, QArrayType, DeltaArrayType, Lambda2ArrayType, LambdaciArrayType> functor(
				jacobianArray, qArray, deltaArray, lambda2Array, lambdaciArray);
			vtkSMPTools::For(0, jacobianArray->GetNumberOfTuples(), functor);
		}
	};
}

void calculateOyf(int xDims, int yDims, int zDims, vtkDataSet* input, vtkFloatArray* oyf_array)
{
	vtkDataArray* vorticityArray = input->GetPointData()->GetArray("vorticity");
	if (vorticityArray == nullptr)
	{
		vorticityArray = input->GetPointData()->GetArray("vorticity-vtk");
	}
	double o_value[3];
	double avg_oy;
	for (int k = 0; k < zDims; k++)
	{
		avg_oy = 0.0;
		// Compute the average oyf value for 1 plane
		for (int j = 0; j < yDims; j++)
		{
			for (int i = 0; i < xDims; i++)
			{
				vtkIdType pointId = i + (xDims * j) + (xDims * yDims * k);
				vorticityArray->GetTuple(pointId, o_value);
				avg_oy += o_value[1]; // Compute Average oy
			}
		}
		avg_oy = avg_oy / (xDims * yDims);

		// Assign a value to oyf_array oyf = oy - avg_oy
		for (int jj = 0; jj < yDims; jj++)
		{
			for (int ii = 0; ii < xDims; ii++)
			{
				vtkIdType pointId = ii + (xDims * jj) + (xDims * yDims * k);
				vorticityArray->GetTuple(pointId, o_value);
				double oy = o_value[1]; // Compute Average oy
				double oyf = oy - avg_oy;
				oyf_array->SetTuple1(pointId, oyf);
			}
		}
	}
}

void calculateVelUf(int xDims, int yDims, int zDims, vtkDataSet* input, 
	vtkFloatArray* vel_uf_array)
{
	vtkDataArray* velocityArray = input->GetPointData()->GetArray("velocity");
	double vel_value[3];
	double avg_u;
	for (int k = 0; k < zDims; k++)
	{
		avg_u = 0.0;
		// Compute the average oyf value for 1 plane
		for (int j = 0; j < yDims; j++)
		{
			for (int i = 0; i < xDims; i++)
			{
				vtkIdType pointId = i + (xDims * j) + (xDims * yDims * k);
				velocityArray->GetTuple(pointId, vel_value);
				avg_u += vel_value[0]; // Compute Average oy
			}
		}
		avg_u = avg_u / (xDims * yDims);

		// Assign a value to oyf_array oyf = oy - avg_oy
		for (int jj = 0; jj < yDims; jj++)
		{
			for (int ii = 0; ii < xDims; ii++)
			{
				vtkIdType pointId = ii + (xDims * jj) + (xDims * yDims * k);
				velocityArray->GetTuple(pointId, vel_value);
				double u = vel_value[0];
				double v = vel_value[0];
				double w = vel_value[0];
				double uF = u - avg_u;
				vel_uf_array->SetTuple3(pointId, uF, v, w);
			}
		}
	}
}

int main(int argc, char* argv[])
{
	vtkSMPTools::SetBackend("STDThread");
	// vtkSMPTools::Initialize(2);
	int num_threads = vtkSMPTools::GetEstimatedNumberOfThreads();
	cout << "Total threads: " << num_threads << endl;

	// std::string datafile = "./fort180_float_bin.vtk";
	std::string datafile = argv[1];
	vtkNew<vtkDataSetReader> dataReader;
	dataReader->SetFileName(datafile.c_str());
	dataReader->Update();

	vtkSmartPointer<vtkStructuredGrid> dataset = vtkStructuredGrid::SafeDownCast(dataReader->GetOutput());
	dataset->Squeeze();

	vtkDataArray* velocity = dataset->GetPointData()->GetArray("velocity");
	if (velocity == nullptr)
	{
		std::cerr << "Velocity array not found." << endl;
		exit(1);
	}

	// 1. Compute the jacobian from the velocity field
	vtkSmartPointer<vtkDataArray> jacobian;
	{
		vtkNew<vtkGradientFilter> gradient;
		gradient->SetInputData(dataset);
		gradient->SetResultArrayName("jacobian-vtk");
		gradient->ComputeDivergenceOn();
		gradient->SetDivergenceArrayName("divergence-vtk");
		gradient->ComputeVorticityOn();
		gradient->SetVorticityArrayName("vorticity-vtk");
		gradient->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, velocity->GetName());
		gradient->Update();

		jacobian = vtkDataArray::SafeDownCast(gradient->GetOutput()->GetPointData()->GetAbstractArray("jacobian-vtk"));
		vtkSmartPointer<vtkDataArray> div = vtkDataArray::SafeDownCast(gradient->GetOutput()->GetPointData()->GetAbstractArray("divergence-vtk"));
		vtkSmartPointer<vtkDataArray> vor = vtkDataArray::SafeDownCast(gradient->GetOutput()->GetPointData()->GetAbstractArray("vorticity-vtk"));
		dataset->GetPointData()->AddArray(jacobian);
		dataset->GetPointData()->AddArray(div);
		dataset->GetPointData()->AddArray(vor);
	}

	// Compute the acceleration field: a = J * v
	vtkSmartPointer<vtkFloatArray> acceleration;
	{
		acceleration = vtkSmartPointer<vtkFloatArray>::New();
		acceleration->SetName("acceleration-vtk");
		acceleration->SetNumberOfComponents(3);
		acceleration->SetNumberOfTuples(jacobian->GetNumberOfTuples());

		MatrixVectorMultiplyWorker worker;

		using Dispatcher = vtkArrayDispatch::Dispatch3ByValueType<vtkArrayDispatch::Reals,
			vtkArrayDispatch::Reals, vtkArrayDispatch::Reals>;

		if (!Dispatcher::Execute(jacobian, velocity, acceleration, worker))
		{
			worker(jacobian.Get(), velocity, acceleration.Get());
		}

		dataset->GetPointData()->AddArray(acceleration);
	}

	vtkSmartPointer<vtkFloatArray> qCriterion;
	vtkSmartPointer<vtkFloatArray> deltaCriterion;
	vtkSmartPointer<vtkFloatArray> lambda2Criterion;
	vtkSmartPointer<vtkFloatArray> lambdaciCriterion;
	{
		qCriterion = vtkSmartPointer<vtkFloatArray>::New();
		deltaCriterion = vtkSmartPointer<vtkFloatArray>::New();
		lambda2Criterion = vtkSmartPointer<vtkFloatArray>::New();
		lambdaciCriterion = vtkSmartPointer<vtkFloatArray>::New();

		qCriterion->SetName("q-vtk");
		deltaCriterion->SetName("delta-vtk");
		lambda2Criterion->SetName("lambda2-vtk");
		lambdaciCriterion->SetName("lambdaci-vtk");

		qCriterion->SetNumberOfComponents(1);
		deltaCriterion->SetNumberOfComponents(1);
		lambda2Criterion->SetNumberOfComponents(1);
		lambdaciCriterion->SetNumberOfComponents(1);

		qCriterion->SetNumberOfTuples(jacobian->GetNumberOfTuples());
		deltaCriterion->SetNumberOfTuples(jacobian->GetNumberOfTuples());
		lambda2Criterion->SetNumberOfTuples(jacobian->GetNumberOfTuples());
		lambdaciCriterion->SetNumberOfTuples(jacobian->GetNumberOfTuples());

		ComputeCriteriaWorker worker;
		worker(jacobian.Get(), qCriterion.Get(), deltaCriterion.Get(), lambda2Criterion.Get(), lambdaciCriterion.Get());

		dataset->GetPointData()->AddArray(qCriterion);
		dataset->GetPointData()->AddArray(deltaCriterion);
		dataset->GetPointData()->AddArray(lambda2Criterion);
		dataset->GetPointData()->AddArray(lambdaciCriterion);
	}

	int dims[3];
	dataset->GetDimensions(dims);

	int xDim = dims[0];
	int yDim = dims[1];
	int zDim = dims[2];

	// vtkNew<vtkIdList> cellsList;
	vtkSmartPointer<vtkFloatArray> vel_uf_array = vtkSmartPointer<vtkFloatArray>::New();
	vel_uf_array->SetName("velocity-uf-vtk");
	vel_uf_array->SetNumberOfComponents(3);
	vel_uf_array->SetNumberOfTuples(jacobian->GetNumberOfTuples());
	calculateVelUf(xDim, yDim, zDim, dataset, vel_uf_array);
	dataset->GetPointData()->AddArray(vel_uf_array);

	// vtkNew<vtkIdList> cellsList;
	vtkSmartPointer<vtkFloatArray> oyf_array = vtkSmartPointer<vtkFloatArray>::New();
	oyf_array->SetName("oyf-vtk");
	oyf_array->SetNumberOfComponents(1);
	oyf_array->SetNumberOfTuples(jacobian->GetNumberOfTuples());
	calculateOyf(xDim, yDim, zDim, dataset, oyf_array);
	dataset->GetPointData()->AddArray(oyf_array);
	
	// Compute & Generate vorticity with oyf
	vtkSmartPointer<vtkFloatArray> vorticity_oyf_array = vtkSmartPointer<vtkFloatArray>::New();
	vorticity_oyf_array->SetName("vorticity-oyf-vtk");
	vorticity_oyf_array->SetNumberOfComponents(3);
	vorticity_oyf_array->SetNumberOfTuples(oyf_array->GetNumberOfTuples());

	vtkDataArray* vorticityArray = dataset->GetPointData()->GetArray("vorticity-vtk");
	float ox, oyf, oz;
	double vorticity_value[3];
	for (int i = 0; i < oyf_array->GetNumberOfTuples(); i++)
	{
		vorticityArray->GetTuple(i, vorticity_value);
		ox = vorticity_value[0];
		oyf = oyf_array->GetTuple1(i);
		oz = vorticity_value[2];
		vorticity_oyf_array->SetTuple3(i, ox, oyf, oz);
	}
	dataset->GetPointData()->AddArray(vorticity_oyf_array);

	vtkNew<vtkGenerateGlobalIds> generator;
	// generator->SetInputData(dataReader->GetOutput());
	generator->SetInputDataObject(dataset);
	generator->Update();

	vtkSmartPointer<vtkDataSet> new_dataset = vtkDataSet::SafeDownCast(generator->GetOutput());
	new_dataset->GetPointData()->SetActiveGlobalIds("");
	new_dataset->GetCellData()->SetActiveGlobalIds("");
	new_dataset->GetPointData()->RemoveArray("vtkGhostType");
	vtkDataArray* arr1 = new_dataset->GetCellData()->GetArray("GlobalCellIds");
	arr1->SetName("DomainCellIds");
	vtkDataArray* arr2 = new_dataset->GetPointData()->GetArray("GlobalPointIds");
	arr2->SetName("DomainPointIds");
	cout << *new_dataset << endl;

	std::string outFileName = argv[2];
	vtkNew<vtkDataSetWriter> writer;
	writer->SetInputData(new_dataset);
	writer->SetFileName(outFileName.c_str());
	writer->SetFileTypeToBinary();
	writer->Write();

	return EXIT_SUCCESS;
}