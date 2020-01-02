/*
        This contains Commandline Interface to the Fractal Phantom
 */
#include "recon_lib.h"

using namespace NDarray;
using namespace std;

int main(int argc, char **argv) {
  // Initialize Phantom
  cout << "Phantom " << endl;
  PHANTOM phantom;
  phantom.read_commandline(argc, argv);
  phantom.init(512, 512, 512, 10);
}
