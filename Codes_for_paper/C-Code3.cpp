#include <Rcpp.h>
using namespace Rcpp;

// Test if the union of "indices" and "k" equals "x"
// Intuition: Test if splitting in the tree "indices" w.r.t. the coordinate "k" is an element of tree "x"
// Input: indices = Indices of a splitable tree, k = splitcoordinate, x = Indices of a tree
// Output: "False" if the union of "indices" and "k" equals "x", otherwise "True"

bool TreeWrong(NumericVector indices, NumericVector k, NumericVector x){
  
  if(is_true(all(in(k,indices)))) {
    
    if(indices.size()==x.size()){
      
      if(is_true(all(in(indices, x)))){
        
        return false;
      }
    }
    
    return true;
  }
  
  if(indices.size()+1==x.size()){
    
    if(is_true(all(in(indices, x)))){
      
      if(is_true(all(in(k,x)))) return false;
    }
  }
  
  return true;
}

// Calculate the optimal split in the current iteration step.
// Input: (X,Y) = data, split_try = number of considered splitpoints in each interval, valiables = coordinates corresponding to currently existing trees,
//        individuals = indices of the data corresponding to the current leaves, leaf_size = minimal leaf size,
//	  split_candidates = coordinates available for splitting in current iteration step
// Output: vector: [0] = minimal achievable sum of squared residuals, [1] = index of the tree, [2] = index of the interval, [3] = coordinate for splitting,
//                 [5] = splitpoint

// [[Rcpp::export]]

NumericVector Calc_Optimal_split2(NumericVector Y, NumericMatrix X, int split_try, List variables, List individuals, List leaf_size, List split_candidates) {
  
  NumericVector R_opt=NumericVector::create(R_PosInf,0,0,0,0);
  
  for(int i_1=0; i_1<variables.size(); ++i_1) {
    
    for(int i_3=0; i_3<split_candidates.size(); ++i_3){
      
      List xxxx=split_candidates(i_3);
      int k= as<int>(xxxx(0))-1;
      
      if(TreeWrong(variables(i_1),xxxx(0),xxxx(1))) continue;
      
      List indiv=individuals(i_1);
      
      for(int i_2=0; i_2<indiv.size(); ++i_2){
        
        NumericVector I=indiv(i_2);
        
        NumericVector Y_1(I.size());
        
        NumericVector Xk_1(I.size());
        
        for(int i_6=0; i_6<I.size(); ++i_6){
          
          Y_1(i_6) = Y(I(i_6)-1);
          
          Xk_1(i_6) =X(I(i_6)-1,k);
        }
        
        NumericVector samplepoints_1=sort_unique(Xk_1);
        
        double leaf_size2 = as<NumericVector>(leaf_size(k))[0];
        
        
        if(samplepoints_1.size()<2*leaf_size2) continue;
        
        NumericVector samplepoints=NumericVector(samplepoints_1.size()+1-2*leaf_size2);
        
        for(int i_5=0; i_5<samplepoints_1.size()+1-2*leaf_size2; ++i_5){
          
          samplepoints(i_5)=samplepoints_1(i_5+leaf_size2);
        }
        
        for(int i_4=0; i_4<split_try; ++i_4){
          
          NumericVector splitpoint=sample(samplepoints,1);
          
          LogicalVector b_1 = Xk_1>=splitpoint(0);
          
          LogicalVector b_2 = Xk_1<splitpoint(0);
          
          NumericVector Y_2=Y_1[b_1];
          
          NumericVector Y_3=Y_1[b_2];
          
          double R=sum(pow(Y_2-mean(Y_2), 2))+sum(pow(Y_3-mean(Y_3) , 2))-sum(pow(Y_2, 2))-sum(pow(Y_3, 2));
          
          if(R_opt(0) <= R) continue;
          
          R_opt(0)=R;
          
          R_opt(1)=i_1+1;
          
          R_opt(2)=i_2+1;
          
          R_opt(3)=k+1;
          
          R_opt(4)=splitpoint(0);
          
        }
      }
    }
  }
  
  return R_opt;
}