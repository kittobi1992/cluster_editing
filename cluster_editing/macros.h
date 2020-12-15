#pragma once

#include <iostream>


template <typename T>
void unused(T&&) {
  // Used to avoid warnings of unused variables
}

// Info, Warning and Error Output Macros
#define GREEN "\033[1;92m"
#define CYAN "\033[1;96m"
#define YELLOW "\033[1;93m"
#define RED "\033[1;91m"
#define BOLD "\033[1m"
#define END "\033[0m"
#define INFO(msg) std::cout << CYAN << "[INFO]" << END << msg << std::endl
#define WARNING(msg) std::cout << YELLOW << "[WARNING]" << END << msg << std::endl
#define ERROR(msg) std::cout << RED << "[ERROR]" << END << msg << std::endl; std::exit(-1)