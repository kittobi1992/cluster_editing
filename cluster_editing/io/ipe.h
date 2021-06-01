/*******************************************************************************
 * This file is part of KaPoCE.
 *
 * KaPoCE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaPoCE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaPoCE.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

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
