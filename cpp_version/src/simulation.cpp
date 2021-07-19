#include <randomPlantedForest.h>
#include <iostream>
#include <cmath>
#include <functional>

double g(std::vector<double> x){
    return sin(sqrt(std::pow(x[0],2)+std::pow(x[1],2)+0.64))/sqrt(std::pow(x[0],2)+std::pow(x[1],2)+0.64);
}

// step-function
double h(std::vector<double> x){
    if(x[0]>0) return 10;
    return -10;
}

// input p: dimension of X_i
// input n: number of data points
void GenerateData(std::vector<double> &Y, std::vector<std::vector<double>> &X, int p=2, int n=100){
	
    std::vector<std::vector<double>> dataPoints;
    //todo: stepSize dependent on n
    for(double i=-10; i<=10; ++i){
        if(p==1){
            dataPoints.push_back({i});
            continue;
        }
        for(double j=-10; j<=10; ++j){
            if(p==2){
                dataPoints.push_back({i, j});
                continue;
            }
            for(double m=-10; m<=10; ++m){
                for(double n=-10; n<=10; ++n){
                    dataPoints.push_back({i, j, m, n});
                }
            }
        }
    }
    std::vector<double> values;
    for(auto point : dataPoints){
        if(p==1) values.push_back(h(point));
        if(p==2) values.push_back(h(point));
        if(p==4) values.push_back(g(point));
    }
    X=dataPoints;
    Y=values;
}


int main() {

    // construct test data
    std::vector<double> Y;
    std::vector<std::vector<double>> X;
    GenerateData(Y, X, 1);

    /*
    for(size_t i = 0; i<Y.size(); ++i){
            for(size_t j = 0; j<X[i].size(); ++j){
                    std::cout << X[i][j] << ',';
            }
            std::cout << Y[i] << std::endl;
    }
    */

    // construct random planted forest
    RandomPlantedForest f1 = RandomPlantedForest(Y, X);
    std::vector<double> res = f1.predict(std::vector<std::vector<double>>{{-10},{-8},{-6},{-4},{-2},{0},{2},{4},{6},{8},{10}});
    //std::vector<double> res = f1.predict(std::vector<std::vector<double>>{{-10,-10},{-8,-8},{-6,-6},{-4,-4},{-2,-2},{0,0},{2,2},{4,4},{6,6},{8,8},{10,10}});
    std::cout << "Predicted Values: "; // << res << std::endl;
    for(auto val: res) std::cout << val << ", ";
    //RandomPlantedForest f2 = RandomPlantedForest(Y, X, 2, 50, 30, std::vector<int> {1,1,1,1}, 10, 0.4, std::vector<std::vector<int>> {{1,2}});

    return(0);
}
