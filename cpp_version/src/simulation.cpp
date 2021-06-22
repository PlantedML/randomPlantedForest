#include <randomPlantedForest.h>
#include <iostream>
#include <cmath>
#include <functional>

double g(std::vector<double> x){
	return sin(sqrt(std::pow(x[0],2)+std::pow(x[1],2)+0.64))/sqrt(std::pow(x[0],2)+std::pow(x[1],2)+0.64);
}

// input p: dimension of X_i
// input n: number of data points
void GenerateData(std::vector<double> &Y, std::vector<std::vector<double>> &X, int p=2, int n=100){
	
    std::vector<std::vector<double>> dataPoints;
    //todo: stepSize dependent on n
    for(double i=-10; i<=10; ++i){
        for(double j=-10; j<=10; ++j){
            for(double m=-10; m<=10; ++m){
                for(double n=-10; n<=10; ++n){
                    dataPoints.push_back({i, j, m, n});
                }
            }
        }
    }
    std::vector<double> values;
    for(auto point : dataPoints){
        values.push_back(g(point));
    }
    X=dataPoints;
    Y=values;
}


int main() {

    // construct test data
    std::vector<double> Y;
    std::vector<std::vector<double>> X;
    GenerateData(Y, X);

    /*
    for(size_t i = 0; i<Y.size(); ++i){
            for(size_t j = 0; j<X[i].size(); ++j){
                    std::cout << X[i][j] << ',';
            }
            std::cout << Y[i] << std::endl;
    }
    */

    // construct random planted forest
    //RandomPlantedForest f1 = RandomPlantedForest(Y, X);
    RandomPlantedForest f2 = RandomPlantedForest(Y, X, 2, 50, 30, std::vector<int> {1,1,1,1}, 10, 0.4, std::vector<std::vector<int>> {{1,2}});

    return(0);
}
