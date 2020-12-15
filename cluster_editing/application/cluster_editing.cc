#include <iostream>

#include "cluster_editing/io/command_line_options.h"

int main(int argc, char* argv[]) {
  cluster_editing::Context context;
  cluster_editing::processCommandLineInput(context, argc, argv);

  std::cout << context << std::endl;
  return 0;
}
