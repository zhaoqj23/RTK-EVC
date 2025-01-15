#pragma once
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include "global.h"

class CombinationGenerator {
public:
	CombinationGenerator(int n, int m) : n(n), m(m) {
		generateCombinations();
	}

	std::vector<std::vector<int>> getCombinations() const {
		return combinations;
	}

private:
	int n;
	int m;
	std::vector<std::vector<int>> combinations;

	void generateCombinations() {
		std::vector<int> combination(m);
		std::vector<int> elements(n);
		std::iota(elements.begin(), elements.end(), 0); // Fill elements with 0, 1, ..., n-1

		generateCombinationsRecursive(elements, combination, 0, 0);
	}

	void generateCombinationsRecursive(const std::vector<int>& elements, std::vector<int>& combination, int start, int depth) {
		if (depth == m) {
			combinations.push_back(combination);
			return;
		}

		for (int i = start; i < n; ++i) {
			combination[depth] = elements[i];
			generateCombinationsRecursive(elements, combination, i + 1, depth + 1);
		}
	}
};
