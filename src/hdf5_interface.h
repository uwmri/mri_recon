#ifndef HDF_INTERFACE
#define HDF_INTERFACE

#include "ArrayTemplates.hpp"
#include <H5Cpp.h>
#include <string>

class HDF5 {
 public:
  // HDF5 File
  H5::H5File file;

  // Constructor
  HDF5(void);
  HDF5(std::string name, const std::string);

  // Allow Adding of Arrays
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<short int, 1> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<short int, 2> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<short int, 3> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<short int, 4> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<short int, 5> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<short int, 6> &A);

  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<float, 1> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<float, 2> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<float, 3> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<float, 4> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<float, 5> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<float, 6> &A);

  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<double, 1> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<double, 2> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<double, 3> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<double, 4> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<double, 5> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<double, 6> &A);

  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<complex<float>, 1> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<complex<float>, 2> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<complex<float>, 3> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<complex<float>, 4> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<complex<float>, 5> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<complex<float>, 6> &A);

  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<complex<double>, 1> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<complex<double>, 2> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<complex<double>, 3> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<complex<double>, 4> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<complex<double>, 5> &A);
  int AddH5Array(std::string GroupName, std::string Name,
                 NDarray::Array<complex<double>, 6> &A);

  // Allow Adding of Arrays
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<float, 1> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<float, 2> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<float, 3> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<float, 4> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<float, 5> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<float, 6> &A);

  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<double, 1> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<double, 2> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<double, 3> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<double, 4> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<double, 5> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<double, 6> &A);

  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<complex<float>, 1> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<complex<float>, 2> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<complex<float>, 3> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<complex<float>, 4> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<complex<float>, 5> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<complex<float>, 6> &A);

  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<complex<double>, 1> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<complex<double>, 2> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<complex<double>, 3> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<complex<double>, 4> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<complex<double>, 5> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<complex<double>, 6> &A);

  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<short, 1> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<short, 2> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<short, 3> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<short, 4> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<short, 5> &A);
  int ReadH5Array(std::string GroupName, std::string Name,
                  NDarray::Array<short, 6> &A);

  int AddH5Scaler(std::string GroupName, std::string Name, int A);
  int AddH5Scaler(std::string GroupName, std::string Name, float A);
  int AddH5Char(std::string GroupName, std::string Name, std::string A);

  int ReadH5Scaler(std::string GroupName, std::string Name, int *A);
  int ReadH5Scaler(std::string GroupName, std::string Name, float *A);
  int ReadH5Char(std::string GroupName, std::string Name, std::string &A);

  static H5::Group open_group(H5::H5File &, std::string GroupName);
  static H5::Group create_group(H5::H5File &, std::string GroupName);

 private:
};

// Template to read a float
template <typename T, int N>
void H5ArrayRead(H5::H5File &file, std::string GroupName, std::string Name,
                 blitz::Array<T, N> &A) {
  // Open the file
  H5::Group group;
  try {
    group = HDF5::open_group(file, GroupName);
  } catch (H5::FileIException &) {
    std::cout << "Can open " << GroupName << " : File Exception" << std::endl;
    throw;
  } catch (H5::GroupIException &) {
    std::cout << "Can open " << GroupName << " : Group Exception" << std::endl;
    throw;
  }

  // Get the dataset
  H5::DataSet dataset;
  try {
    dataset = group.openDataSet(Name);
  } catch (...) {
    std::cout << "Can't open Group " << GroupName << " / " << Name << std::endl;
    throw;
  }
  // Get the dataspace
  H5::DataSpace dataspace = dataset.getSpace();

  // Get the number of dimensions
  int rank = dataspace.getSimpleExtentNdims();
  if (rank != N) {
    std::cout << "Expected Rank " << N << " but opened a rank " << rank
              << " file " << std::endl
              << std::flush;
  }

  // Create DataSet
  hsize_t dimsf[6];
  int ndims = dataspace.getSimpleExtentDims(dimsf, NULL);

  unsigned int total = dataspace.getSelectNpoints();
  if (total != A.numElements()) {
    std::cout << "Dataspace has " << total << " points while array has "
              << A.numElements() << std::endl;
    std::cout << "Resizing input array" << std::endl;
    NDarray::TinyVector<int, N> dims;
    for (int i = 0; i < N; i++) {
      dims(i) = dimsf[ndims - i - 1];
    }
    A.setStorage(NDarray::ColumnMajorArray<N>());
    A.resize(dims);
  }

  H5::DataType data_type = H5::DataType(dataset.getDataType());

  // Read the data
  dataset.read(A.data(), data_type);
}

// Template to write a complex<double>
template <int N>
void H5ArrayWrite(H5::H5File &file, std::string GroupName, std::string Name,
                  blitz::Array<complex<double>, N> &A) {
  H5::Exception::dontPrint();
  H5::Group group = HDF5::create_group(file, GroupName);

  // Create DataSet
  hsize_t dimsf[N];  // dataset dimensions
  for (int i = 0; i < N; i++) {
    dimsf[i] = A.length(N - i - 1);
  }
  H5::DataSpace dataspace(N, dimsf);

  /*DataType Needed for Complex<float>*/
  H5::CompType datatype(sizeof(complex<float>));
  datatype.insertMember("real", 0, H5::PredType::NATIVE_DOUBLE);
  datatype.insertMember("imag", sizeof(float), H5::PredType::NATIVE_DOUBLE);

  H5::DataSet dataset(group.createDataSet(Name, datatype, dataspace));
  dataset.write(A.data(), datatype, dataspace);
}

// Template to write a complex<float>
template <int N>
void H5ArrayWrite(H5::H5File &file, std::string GroupName, std::string Name,
                  blitz::Array<complex<float>, N> &A) {
  H5::Exception::dontPrint();
  H5::Group group = HDF5::create_group(file, GroupName);

  // Create DataSet
  hsize_t dimsf[N];  // dataset dimensions
  for (int i = 0; i < N; i++) {
    dimsf[i] = A.length(N - i - 1);
  }
  H5::DataSpace dataspace(N, dimsf);

  /*DataType Needed for Complex<float>*/
  H5::CompType datatype(sizeof(complex<float>));
  datatype.insertMember("real", 0, H5::PredType::NATIVE_FLOAT);
  datatype.insertMember("imag", sizeof(float), H5::PredType::NATIVE_FLOAT);

  H5::DataSet dataset(group.createDataSet(Name, datatype, dataspace));
  dataset.write(A.data(), datatype, dataspace);
}

// Template to write a double
template <int N>
void H5ArrayWrite(H5::H5File &file, std::string GroupName, std::string Name,
                  blitz::Array<double, N> &A) {
  H5::Exception::dontPrint();
  H5::Group group = HDF5::create_group(file, GroupName);

  // Create DataSet
  hsize_t dimsf[N];  // dataset dimensions
  for (int i = 0; i < N; i++) {
    dimsf[i] = A.length(N - i - 1);
  }
  H5::DataSpace dataspace(N, dimsf);

  H5::DataSet dataset(
      group.createDataSet(Name, H5::PredType::NATIVE_DOUBLE, dataspace));
  dataset.write(A.data(), H5::PredType::NATIVE_DOUBLE, dataspace);
}

// Template to write a float
template <int N>
void H5ArrayWrite(H5::H5File &file, std::string GroupName, std::string Name,
                  NDarray::Array<float, N> &A) {
  H5::Exception::dontPrint();
  H5::Group group = HDF5::create_group(file, GroupName);

  // Create DataSet
  hsize_t dimsf[N];  // dataset dimensions
  for (int i = 0; i < N; i++) {
    dimsf[i] = A.length(N - i - 1);
  }
  H5::DataSpace dataspace(N, dimsf);

  H5::DataSet dataset(
      group.createDataSet(Name, H5::PredType::NATIVE_FLOAT, dataspace));
  dataset.write(A.data(), H5::PredType::NATIVE_FLOAT, dataspace);
}

// Template to write a short int
template <int N>
void H5ArrayWrite(H5::H5File &file, std::string GroupName, std::string Name,
                  NDarray::Array<short, N> &A) {
  H5::Exception::dontPrint();
  H5::Group group = HDF5::create_group(file, GroupName);

  // Create DataSet
  hsize_t dimsf[N];  // dataset dimensions
  for (int i = 0; i < N; i++) {
    dimsf[i] = A.length(N - i - 1);
  }
  H5::DataSpace dataspace(N, dimsf);

  H5::DataSet dataset(
      group.createDataSet(Name, H5::PredType::NATIVE_SHORT, dataspace));
  dataset.write(A.data(), H5::PredType::NATIVE_SHORT, dataspace);
}

#endif
