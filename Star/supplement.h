#include <boundaries.h>
#include <fstream>

#ifndef SUPPLEMENT_H
#define SUPPLEMENT_H

Domain_Partitioning read_partitioning(const std::string &file_name);
std::ostream &operator<<(std::ostream &os, const Domain_Partitioning &part);

#endif