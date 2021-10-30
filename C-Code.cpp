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

NumericVector Calc_Optimal_split2(NumericVector Y, NumericVector W, NumericMatrix X, int split_try, List variables, List individuals, List leaf_size, List split_candidates, StringVector loss, std::vector<int>  categorical_columns, int max_categorical, double delta) {

  String loss2(loss[0]);

  NumericVector R_opt (max_categorical+6, R_PosInf);


  for(int i_1=0; i_1<variables.size(); ++i_1) {

    for(int i_3=0; i_3<split_candidates.size(); ++i_3){

      List xxxx=split_candidates(i_3);
      int k= as<int>(xxxx(0))-1;


      if(TreeWrong(variables(i_1),xxxx(0),xxxx(1))) continue;

      List indiv=individuals(i_1);

      for(int i_2=0; i_2<indiv.size(); ++i_2){

        NumericVector I=indiv(i_2);

        NumericVector Y_1(I.size());

        NumericVector W_1(I.size());


        NumericVector Xk_1(I.size());

        for(int i_6=0; i_6<I.size(); ++i_6){

          Y_1(i_6) = Y(I(i_6)-1);

          W_1(i_6) = W(I(i_6)-1);

          Xk_1(i_6) =X(I(i_6)-1,k);
        }

        NumericVector samplepoints_1=sort_unique(Xk_1);         /// Number of unique values leaf

        double leaf_size2 = as<NumericVector>(leaf_size(k))[0];  /// minimal Number of individuals that should be in a leaf for variable k


        if(samplepoints_1.size()<2*leaf_size2) continue;





        for(int i_4=0; i_4<split_try; ++i_4){

          LogicalVector b_1 (Xk_1.size());
          LogicalVector b_2 (Xk_1.size());
          NumericVector splitpoint (1);

          bool is_categorical = std::find(std::begin(categorical_columns), std::end(categorical_columns), k+1) != std::end(categorical_columns);

          if(is_categorical) {
            Rcpp::IntegerVector group_sizes = Rcpp::seq(1, samplepoints_1.size());
            double group_size = Rcpp::sample(group_sizes,1)(0);
            NumericVector group1 = Rcpp::sample(samplepoints_1, group_size);



            for(int i_x=0; i_x < Xk_1.size(); ++i_x){
              b_1(i_x)  = std::find(std::begin(group1), std::end(group1), Xk_1(i_x)) != std::end(group1);
            };
            b_2 = !b_1;

          } else {

            NumericVector samplepoints=NumericVector(samplepoints_1.size()+1-2*leaf_size2);

            for(int i_5=0; i_5<samplepoints_1.size()+1-2*leaf_size2; ++i_5){
              samplepoints(i_5)=samplepoints_1(i_5+leaf_size2);
            };

          splitpoint=sample(samplepoints,1);

          b_1 = Xk_1>=splitpoint(0);

          b_2 = Xk_1<splitpoint(0);
          };

          NumericVector v (I.size());

         // if (loss2=="exponential"){
         //   LogicalVector bb_1 = W_1>v;
        //    LogicalVector bb_2 = W_1>v;

        //    b_1 = b_1&&bb_1;
        //    b_2 = b_2&&bb_2;
      //    }

          NumericVector Y_2=Y_1[b_1];

          NumericVector Y_3=Y_1[b_2];

          NumericVector W_2=W_1[b_1];

          NumericVector W_3=W_1[b_2];

          double R;



          if ( loss2=="L2") {
            R=sum(pow(Y_2-mean(Y_2), 2))+sum(pow(Y_3-mean(Y_3) , 2))-sum(pow(Y_2, 2))-sum(pow(Y_3, 2));
          } else if ( loss2=="L1"){ R=sum(abs(Y_2-mean(Y_2)))+sum(abs(Y_3-mean(Y_3)))-sum(abs(Y_2))-sum(abs(Y_3));
	  } else if ( loss2=="median"){ R=sum(abs(Y_2-median(Y_2)))+sum(abs(Y_3-median(Y_3)))-sum(abs(Y_2))-sum(abs(Y_3));
          } else if ( loss2=="exponential"){



            double R21;
            double R22;
            double R31;
            double R32;

              R21 = sum(((Y_2+1)/2) *(W_2/sum(W_2)));
              R31 = sum(((Y_3+1)/2) *(W_3/sum(W_3)));
              R21 = std::min(1-delta,std::max(delta,R21));
              R22 = 1-R21;
              R31 = std::min(1-delta,std::max(delta,R31));
              R32 = 1-R31;


              R= sum(W_2*exp(-0.5*Y_2*log( R21/R22 ))) + sum(W_3*exp(-0.5*Y_3*log( R31/R32 ))) - sum(W_2)-sum(W_3);

             if(R==0) {
                R = R_PosInf;
            };
            if(std::isnan(R)==TRUE){
              R = R_PosInf;
            };
       //       if((sum(W_2*exp(-0.5*Y_2*log( mean(((Y_2+1)/2) *(W_2))/(mean(((Y_2-1)/-2)*(W_2))) )))==0)&&(sum(exp(-0.5*Y_2*log( mean(((Y_2+1)/2) )/(mean(((Y_2-1)/-2))) )))!=0)){
      //          R = R_PosInf;
      //       };
       //       if((sum(W_3*exp(-0.5*Y_3*log( mean(((Y_3+1)/2) *(W_3))/(mean(((Y_3-1)/-2)*(W_3))) )))==0)&&(sum(exp(-0.5*Y_3*log( mean(((Y_3+1)/2) )/(mean(((Y_3-1)/-2))) )))!=0)){
       //        R = R_PosInf;
       //      };


          }else {

            NumericVector P_1 = 1/(1+exp(-W_1));

            double R_old = sum(Y_1*log(P_1) + (1-Y_1)*log(1-P_1));

            double R2 = mean(Y_2);
            double R3 = mean(Y_3);

            R2 = std::min(1-delta,std::max(delta, R2));
            R3 = std::min(1-delta,std::max(delta, R3));

            W_2 = W_2 +  log(R2/(1-R2)) - mean(W_2);
            W_3 = W_3 +  log(R3/(1-R3)) - mean(W_3);


            NumericVector P_2 = 1/(1+exp(-W_2));
            NumericVector P_3 = 1/(1+exp(-W_3));

            R = -sum(Y_2*log(P_2)+ (1-Y_2)*log(1-P_2)) - sum(Y_3*log(P_3)+ (1-Y_3)*log(1-P_3)) + R_old;

            if(std::isnan(R)==TRUE){
              R = R_PosInf;
            };


          };


          if(R_opt(0) <= R) continue;

          R_opt(0)=R;

          R_opt(1)=i_1+1;

          R_opt(2)=i_2+1;

          R_opt(3)=k+1;

          R_opt(5)=1;

          if(!is_categorical){ R_opt(4)=splitpoint(0);
            for (int i=0; i< max_categorical-2; ++i){
              R_opt(i+6) = R_PosInf;
            };
                        };

          if(is_categorical){
            R_opt(5)=0;
            NumericVector temp = Xk_1[b_1];
            NumericVector splitpoints_2 = sort_unique(temp);
            for (int i=0; i< max_categorical-2 ; ++i){
              if (i < splitpoints_2.size()){
              R_opt(i+6) = splitpoints_2(i);} else {
                R_opt(i+6) = R_PosInf;
              };

            };
          };



        }
      }
    }
  }

  return R_opt;
}
