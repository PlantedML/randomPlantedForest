#include <random_planted_forest.h>
#include <iostream>
#include <math.h>
#include <functional>

double g(std::vector<double> x){
	return sin(sqrt(std::pow(x[0],2)+std::pow(x[1],2)+0.64))/sqrt(std::pow(x[0],2)+std::pow(x[1],2)+0.64);
}

// input p: dimension of X_i
// input n: number of data points
void GenerateData(std::vector<double> &Y, std::vector<std::vector<double>> &X, int p=2, int n=100){
	std::vector<std::vector<double>> dataPoints;

	//todo: stepSize dependent on n 
	for(int i=-10; i<10; ++i){
		std::vector<double> vec;
		for(int j=-10; j<10; ++j){
			vec={i,j};
		}
		dataPoints.push_back(vec);
	}
	std::vector<double> values;
	for(auto point : dataPoints){
		values.push_back(g(point));
	}
	X=dataPoints;
	Y=values;
}


int main() {

	std::vector<double> Y;
	std::vector<std::vector<double>> X;
	GenerateData(Y, X);
	for(auto el:Y){
		std::cout << el << std::endl;
	}
	RandomPlantedForest f1 = RandomPlantedForest(Y, X);

    return(0);
}
