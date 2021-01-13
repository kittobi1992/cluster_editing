#include "cluster_editing/io/ipe.h"

#include <iomanip>
#include <ios>
#include <string>

IpeFile::IpeFile(const std::string& filename) : _file(filename) {
  _file << std::fixed << std::setprecision(1);
  file_start();
  page_start();
}

IpeFile::~IpeFile() {
  page_end();
  file_end();
}

void IpeFile::new_page() {
  page_end();
  page_start();
}

void IpeFile::label(const std::string& label, double x, double y,
                    const std::string& color) {
  _file << R"(<text layer="alpha" transformations="translations" pos=")" << x
         << " " << y << "\" stroke=\"" << color << R"(" type="label" )"
         << R"(halign="center" size="normal" valign="center">)" << label
         << "</text>\n";
}

void IpeFile::line(double x1, double y1, double x2, double y2,
                   const std::string& color, double pen, bool transparent) {
  _file << "<path stroke = \"" << color << "\" pen = \"" << pen << "\""
        << (transparent ? "stroke-opacity=\"transparent\"" : "") << ">\n";
  _file << x1 << " " << y1 << " m\n";
  _file << x2 << " " << y2 << " l\n";
  _file << "</path>\n";
}

void IpeFile::point(double x, double y, const std::string& color) {
  _file << "<use name=\"mark/disk(sx)\" pos=\"" << x << " " << y
        << R"(" size="normal" stroke=")" << color << "\"/>\n";
}

void IpeFile::disk(double x, double y, double radius,
                   const std::string& color) {
  _file << "<path fill=\"" << color << "\" opacity=\"transparent\">\n"
        << radius << " 0 0 " << radius << " " << x << " " << y << " e\n"
        << "</path>\n";
}

void IpeFile::file_start() {
  _file << "<?xml version=\"1.0\"?>\n"
        << "<!DOCTYPE ipe SYSTEM \"ipe.dtd\">\n"
        << "<ipe version=\"70206\" creator=\"Ipe 7.2.6\">\n";
  _file << "<ipestyle name=\"dummy\">\n";
  _file << "<symbol name=\"mark/disk(sx)\" transformations=\"translations\">\n"
        << "<path fill=\"sym-stroke\">\n"
        << "0.6 0 0 0.6 0 0 e\n"
        << "</path>\n"
        << "</symbol>\n";
  _file << "<opacity name=\"transparent\" value=\"0.4\"/>\n";
  _file << "<textstretch name=\"normal\" value=\"0.2\"/>\n";
  _file << "<layout paper=\"2000 2000\" origin=\"0 0\" frame=\"2000 2000\"/>\n";
  _file << "</ipestyle>\n";
}

void IpeFile::file_end() { _file << "</ipe>\n"; }

void IpeFile::page_start() {
  _file << "<page>\n";
  _file << "<layer name=\"alpha\"/>\n";
  _file << "<view layers=\"alpha\" active=\"alpha\"/>\n";
}

void IpeFile::page_end() { _file << "</page>\n"; }
