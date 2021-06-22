#ifndef DECISION_TREE_H
#define DECISION_TREE_H

#include <vector>
#include <utility>
#include <set>

typedef std::pair<double, double> Interval;

struct Leaf{
    std::set<int> individuals;          // considered samples for each leaf
    double value;                       // residual
    std::vector<Interval> intervals;    // min/max for each feature of the interval
};

class DecisionTree {

    friend class RandomPlantedForest;

    public:


    private:
        std::set<int> split_dims;           // dimensions of the performed splits
        std::vector<Leaf> leaves;           // leaves of tree containing intervals and approximating value
        // idea: save intervals as interval-tree with nodes and corresponding values
};


#endif // DECISION_TREE_H
