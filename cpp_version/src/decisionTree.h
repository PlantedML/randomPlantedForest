#ifndef DECISION_TREE_H
#define DECISION_TREE_H

#include <vector>
#include <utility>
#include <set>

struct Split {
	double min_sum;			//
	int tree_index;			//
	int interval_index;		//
	int split_coordinate;	//
	double split_point;		//
};

class DecisionTree {
	
	public:
		
	private:
		std::set<int> split_dims; 			// dimensions of the performed splits
		std::vector<double> values;			// residuals of the leaves
		std::vector<std::pair<double, double>> intervals;	// min/max of interval for each leaf
		Split calcOptimalSplit();			// determine optimal split
		
};


#endif // DECISION_TREE_H
