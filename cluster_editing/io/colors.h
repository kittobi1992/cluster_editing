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
