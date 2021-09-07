#ifndef _OUTPUT_H
#define _OUTPUT_H

#include <vector>
#include <fstream>
#include <memory>
#include "World.h"
#include "Species.h"

namespace Output {
	void fields(World &world, std::vector<std::unique_ptr<Species>> &species);
	void screenOutput(World &world, std::vector<std::unique_ptr<Species>> &species);
	void diagOutput(World &world, std::vector<std::unique_ptr<Species>> &species);
	void particles(World &world, std::vector<std::unique_ptr<Species>> &species, int num_parts);
}

#endif
