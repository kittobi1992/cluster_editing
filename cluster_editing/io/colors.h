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

#include <sstream>
#include <string>

struct Color {
  double r;
  double g;
  double b;

  operator std::string() const {
    std::stringstream result;
    result << r << " " << g << " " << b;
    return result.str();
  }
};

Color rgb(double R, double G, double B);
Color rgb255(unsigned R, unsigned G, unsigned B);
Color hsv(unsigned H, double S, double V);
