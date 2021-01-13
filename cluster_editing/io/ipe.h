#pragma once

#include <fstream>
#include <string>

#include "cluster_editing/io/colors.h"

class IpeFile {
 public:
  explicit IpeFile(const std::string& filename);
  ~IpeFile();

  void new_page();

  void label(const std::string& label, double x, double y,
             const std::string& color = "black");

  void line(double x1, double y1, double x2, double y2,
            const std::string& color = "black", double pen = 0.4,
            bool transparent = false);

  void point(double x, double y, const std::string& color = "black");

  void disk(double x, double y, double radius,
            const std::string& color = "black");

 private:
  std::ofstream _file;
  
  void file_start();
  void file_end();
  void page_start();
  void page_end();
};
