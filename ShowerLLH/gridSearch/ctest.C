//#include <fstream>
//#include <iomanip>
//#include <sstream>
//#include <string>
#include <vector>

using namespace std;

//void testFunc(vector<double> grid, 
//    vector< vector< vector< vector< vector<double> > > > > llhTable, 
//    double Zbin, vector<double> Sbin, vector<double> Cbin, 
//    double theta, double phi, vector< vector<double> > tankxyz, 
//    vector<double> Dbins) {

//void testFunc(vector<double> grid, 
//    vector<double> llhTable) {

//  int nG = grid.size();
//  int nE = llhTable.size();

//  cout << nG << endl;
//  cout << nE << endl;
//}

char const* testFunc() {
  return "Hello, world";
}

#include <boost/python.hpp>

BOOST_PYTHON_MODULE(ctest) {
  using namespace boost::python;
  def("testFunc", testFunc);
}
