#include "hdf5_interface.h"

using namespace std;
using namespace NDarray;

HDF5::HDF5(string FileName, const string operation) {
  cout << "HDF5 - Opening file operation = " << operation << endl;

  // OPen the file
  if (operation.compare("r") == 0) {
    file = H5::H5File(FileName, H5F_ACC_RDONLY);
  } else if ((operation.compare("rw") == 0) || (operation.compare("rw") == 0)) {
    file = H5::H5File(FileName, H5F_ACC_RDWR);
  } else {
    file = H5::H5File(FileName, H5F_ACC_TRUNC);
  }
}

HDF5::HDF5(void) {}

H5::Group HDF5::create_group(H5::H5File &act_file, string GroupName) {
  // Create Group
  H5::Group group;
  try {
    // cout << "Trying to Open: " << GroupName << endl;
    group = H5::Group(act_file.openGroup(GroupName));
  } catch (H5::FileIException &error) {
    // cout << "Group does not exist-Creating" << endl;
    string x(GroupName);
    x.insert(0, "/");
    cout << x << endl;
    group = H5::Group(act_file.createGroup(x));
  }

  return (group);
}

H5::Group HDF5::open_group(H5::H5File &act_file, string GroupName) {
  // Open Group
  H5::Group group;
  try {
    group = H5::Group(act_file.openGroup(GroupName));
  } catch (H5::FileIException &error) {
    // cout << "Group does not exist-Creating" << endl;
    string x(GroupName);
    x.insert(0, "/");
    cout << x << endl;
    cout << "HDF5 Read : Group " << GroupName << " does not exist" << endl;
  }
  return (group);
}

int HDF5::ReadH5Scaler(string GroupName, string Name, int *out) {
  H5::Exception::dontPrint();

  // Open Group
  H5::Group group = open_group(file, GroupName);

  // Now open attribute
  H5::Attribute att;
  try {
    att = group.openAttribute(Name);
  } catch (H5::FileIException &error) {
    cout << "HDF5 Read :  " << GroupName << " attribute " << Name
         << " does not exist" << endl;
  }

  // Get the data type
  H5::DataType data_type = H5::DataType(att.getDataType());
  H5::DataType int_type(H5::PredType::STD_I32LE);

  if (int_type == data_type) {
    att.read(data_type, out);
  } else {
    cout << "Expected to read little endian int32 but did not" << endl;
    exit(1);
  }

  return (0);
}

int HDF5::AddH5Scaler(string GroupName, string Name, int A) {
  H5::Exception::dontPrint();

  // Create Group
  H5::Group group = create_group(file, GroupName);

  H5::IntType int_type(H5::PredType::STD_I32LE);
  H5::DataSpace att_space(H5S_SCALAR);
  H5::Attribute att = group.createAttribute(Name, int_type, att_space);
  att.write(int_type, &A);

  return (0);
}

int HDF5::ReadH5Scaler(string GroupName, string Name, float *out) {
  H5::Exception::dontPrint();

  // Open Group
  H5::Group group = open_group(file, GroupName);

  // Now open attribute
  H5::Attribute att;
  try {
    att = group.openAttribute(Name);
  } catch (H5::FileIException &error) {
    cout << "HDF5 Read :  " << GroupName << " attribute " << Name
         << " does not exist" << endl;
  }

  // Get the data type
  H5::DataType data_type = H5::DataType(att.getDataType());
  H5::DataType int_type(H5::PredType::NATIVE_FLOAT);

  if (int_type == data_type) {
    att.read(data_type, out);
  } else {
    cout << "Expected to read little endian float32 but did not" << endl;
    exit(1);
  }

  return (0);
}

int HDF5::AddH5Scaler(string GroupName, string Name, float A) {
  H5::Exception::dontPrint();

  // Create Group
  H5::Group group = create_group(file, GroupName);

  H5::FloatType float_type(H5::PredType::NATIVE_FLOAT);
  H5::DataSpace att_space(H5S_SCALAR);
  H5::Attribute att = group.createAttribute(Name, float_type, att_space);
  att.write(float_type, &A);

  return (0);
}

int HDF5::ReadH5Char(string GroupName, string Name, string &out) {
  H5::Exception::dontPrint();

  // Open Group
  H5::Group group = open_group(file, GroupName);

  // Now open attribute
  H5::Attribute att;
  try {
    att = group.openAttribute(Name);
  } catch (H5::FileIException &error) {
    cout << "HDF5 Read :  " << GroupName << " attribute " << Name
         << " does not exist" << endl;
  }

  // Get the data type
  H5std_string strreadbuf("");
  H5::StrType str_type(H5::PredType::C_S1, 512);
  att.read(str_type, strreadbuf);

  out = strreadbuf;

  return (0);
}

int HDF5::AddH5Char(string GroupName, string Name, string S) {
  H5::Exception::dontPrint();

  // Create Group
  H5::Group group = create_group(file, GroupName);

  H5::StrType str_type(H5::PredType::C_S1, 512);
  H5::DataSpace att_space(H5S_SCALAR);

  H5::Attribute att = group.createAttribute(Name, str_type, att_space);
  att.write(str_type, S);

  return (0);
}

int HDF5::AddH5Array(string GroupName, string Name, Array<short, 1> &A) {
  H5ArrayWrite<1>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<short, 2> &A) {
  H5ArrayWrite<2>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<short, 3> &A) {
  H5ArrayWrite<3>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<short, 4> &A) {
  H5ArrayWrite<4>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<short, 5> &A) {
  H5ArrayWrite<5>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<short, 6> &A) {
  H5ArrayWrite<6>(file, GroupName, Name, A);
  return (0);
}

int HDF5::AddH5Array(string GroupName, string Name, Array<float, 1> &A) {
  H5ArrayWrite<1>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<float, 2> &A) {
  H5ArrayWrite<2>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<float, 3> &A) {
  H5ArrayWrite<3>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<float, 4> &A) {
  H5ArrayWrite<4>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<float, 5> &A) {
  H5ArrayWrite<5>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<float, 6> &A) {
  H5ArrayWrite<6>(file, GroupName, Name, A);
  return (0);
}

int HDF5::AddH5Array(string GroupName, string Name, Array<u_int64_t, 1> &A) {
  H5ArrayWrite<1>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<u_int64_t, 2> &A) {
  H5ArrayWrite<2>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<u_int64_t, 3> &A) {
  H5ArrayWrite<3>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<u_int64_t, 4> &A) {
  H5ArrayWrite<4>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<u_int64_t, 5> &A) {
  H5ArrayWrite<5>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<u_int64_t, 6> &A) {
  H5ArrayWrite<6>(file, GroupName, Name, A);
  return (0);
}

int HDF5::AddH5Array(string GroupName, string Name, Array<double, 1> &A) {
  H5ArrayWrite<1>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<double, 2> &A) {
  H5ArrayWrite<2>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<double, 3> &A) {
  H5ArrayWrite<3>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<double, 4> &A) {
  H5ArrayWrite<4>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<double, 5> &A) {
  H5ArrayWrite<5>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name, Array<double, 6> &A) {
  H5ArrayWrite<6>(file, GroupName, Name, A);
  return (0);
}

int HDF5::AddH5Array(string GroupName, string Name,
                     Array<complex<float>, 1> &A) {
  H5ArrayWrite<1>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name,
                     Array<complex<float>, 2> &A) {
  H5ArrayWrite<2>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name,
                     Array<complex<float>, 3> &A) {
  H5ArrayWrite<3>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name,
                     Array<complex<float>, 4> &A) {
  H5ArrayWrite<4>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name,
                     Array<complex<float>, 5> &A) {
  H5ArrayWrite<5>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name,
                     Array<complex<float>, 6> &A) {
  H5ArrayWrite<6>(file, GroupName, Name, A);
  return (0);
}

int HDF5::AddH5Array(string GroupName, string Name,
                     NDarray::Array<complex<double>, 1> &A) {
  H5ArrayWrite<1>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name,
                     NDarray::Array<complex<double>, 2> &A) {
  H5ArrayWrite<2>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name,
                     NDarray::Array<complex<double>, 3> &A) {
  H5ArrayWrite<3>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name,
                     NDarray::Array<complex<double>, 4> &A) {
  H5ArrayWrite<4>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name,
                     NDarray::Array<complex<double>, 5> &A) {
  H5ArrayWrite<5>(file, GroupName, Name, A);
  return (0);
}
int HDF5::AddH5Array(string GroupName, string Name,
                     NDarray::Array<complex<double>, 6> &A) {
  H5ArrayWrite<6>(file, GroupName, Name, A);
  return (0);
}

int HDF5::ReadH5Array(string GroupName, string Name, Array<float, 1> &A) {
  H5ArrayRead<float, 1>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name, Array<float, 2> &A) {
  H5ArrayRead<float, 2>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name, Array<float, 3> &A) {
  H5ArrayRead<float, 3>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name, Array<float, 4> &A) {
  H5ArrayRead<float, 4>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name, Array<float, 5> &A) {
  H5ArrayRead<float, 5>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name, Array<float, 6> &A) {
  H5ArrayRead<float, 6>(file, GroupName, Name, A);
  return (0);
}

int HDF5::ReadH5Array(string GroupName, string Name, Array<u_int64_t, 1> &A) {
  H5ArrayRead<u_int64_t, 1>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name, Array<u_int64_t, 2> &A) {
  H5ArrayRead<u_int64_t, 2>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name, Array<u_int64_t, 3> &A) {
  H5ArrayRead<u_int64_t, 3>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name, Array<u_int64_t, 4> &A) {
  H5ArrayRead<u_int64_t, 4>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name, Array<u_int64_t, 5> &A) {
  H5ArrayRead<u_int64_t, 5>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name, Array<u_int64_t, 6> &A) {
  H5ArrayRead<u_int64_t, 6>(file, GroupName, Name, A);
  return (0);
}

int HDF5::ReadH5Array(string GroupName, string Name, Array<double, 1> &A) {
  H5ArrayRead<double, 1>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name, Array<double, 2> &A) {
  H5ArrayRead<double, 2>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name, Array<double, 3> &A) {
  H5ArrayRead<double, 3>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name, Array<double, 4> &A) {
  H5ArrayRead<double, 4>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name, Array<double, 5> &A) {
  H5ArrayRead<double, 5>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name, Array<double, 6> &A) {
  H5ArrayRead<double, 6>(file, GroupName, Name, A);
  return (0);
}

int HDF5::ReadH5Array(string GroupName, string Name,
                      Array<complex<float>, 1> &A) {
  H5ArrayRead<complex<float>, 1>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name,
                      Array<complex<float>, 2> &A) {
  H5ArrayRead<complex<float>, 2>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name,
                      Array<complex<float>, 3> &A) {
  H5ArrayRead<complex<float>, 3>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name,
                      Array<complex<float>, 4> &A) {
  H5ArrayRead<complex<float>, 4>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name,
                      Array<complex<float>, 5> &A) {
  H5ArrayRead<complex<float>, 5>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name,
                      Array<complex<float>, 6> &A) {
  H5ArrayRead<complex<float>, 6>(file, GroupName, Name, A);
  return (0);
}

int HDF5::ReadH5Array(string GroupName, string Name,
                      Array<complex<double>, 1> &A) {
  H5ArrayRead<complex<double>, 1>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name,
                      Array<complex<double>, 2> &A) {
  H5ArrayRead<complex<double>, 2>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name,
                      Array<complex<double>, 3> &A) {
  H5ArrayRead<complex<double>, 3>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name,
                      Array<complex<double>, 4> &A) {
  H5ArrayRead<complex<double>, 4>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name,
                      Array<complex<double>, 5> &A) {
  H5ArrayRead<complex<double>, 5>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name,
                      Array<complex<double>, 6> &A) {
  H5ArrayRead<complex<double>, 6>(file, GroupName, Name, A);
  return (0);
}

int HDF5::ReadH5Array(string GroupName, string Name,
                      Array<short, 1> &A) {
  H5ArrayRead<short, 1>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name,
                      Array<short, 2> &A) {
  H5ArrayRead<short, 2>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name,
                      Array<short, 3> &A) {
  H5ArrayRead<short, 3>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name,
                      Array<short, 4> &A) {
  H5ArrayRead<short, 4>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name,
                      Array<short, 5> &A) {
  H5ArrayRead<short, 5>(file, GroupName, Name, A);
  return (0);
}
int HDF5::ReadH5Array(string GroupName, string Name,
                      Array<short, 6> &A) {
  H5ArrayRead<short, 6>(file, GroupName, Name, A);
  return (0);
}
