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

#include "cluster_editing/io/colors.h"

#include <cmath>

Color rgb(double R, double G, double B) { return Color{R, G, B}; }

Color hsv(unsigned H, double S, double V) {
  H = H % 360;
  double C = V * S;
  double X = C * (1 - std::abs(std::fmod((H / 60.0), 2) - 1));
  double m = V - C;
  double R = 0, G = 0, B = 0;
  if (H < 60) {
    R = C;
    G = X;
    B = 0;
  } else if (60 <= H && H < 120) {
    R = X;
    G = C;
    B = 0;
  } else if (120 <= H && H < 180) {
    R = 0;
    G = C;
    B = X;
  } else if (180 <= H && H < 240) {
    R = 0;
    G = X;
    B = C;
  } else if (240 <= H && H < 300) {
    R = X;
    G = 0;
    B = C;
  } else if (300 <= H && H < 360) {
    R = C;
    G = 0;
    B = X;
  }
  return rgb(R + m, G + m, B + m);
}
