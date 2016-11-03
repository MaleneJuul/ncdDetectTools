#include <RcppArmadillo.h>
#include <armadillo>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]


mat convolution_body(mat dataset, int threshold, vec integers) {
    
    // construct needed vectors and integers
    vec p_out = zeros<vec> (threshold + 1) ;
    vec p_tmp = zeros<vec> (threshold + 1) ;
    
    vec y = dataset.col(0) ;
    vec p = dataset.col(1) ;
    vec x = dataset.col(2) ;
    
    int max_x = x.max() ;
    int min_x = x.min() ;
    
    // set probabilities for the first x-value
    uvec indices = find(x == min_x) ;
    
    vec yi = y.elem(indices) ;
    vec pi = p.elem(indices) ;
    
    for (int j = 0; j < yi.n_elem; j = j + 1 ) {
        if (yi[j] <= threshold) {
            
            // take existing probability for current y, and add new value
            p_out[yi[j]] = p_out[yi[j]] + pi[j];
            
        } else {
            
            // if value greater than threshold, then place in last slot of output vector
            p_out[threshold] = p_out[threshold] + pi[j];
            
        }
    }
    
    
    uvec non_zero_ps = find(p_out != 0) ;
    vec non_zero_ys = integers(non_zero_ps) ;
    
    // go through the remaining x-values
    for (int i = (min_x + 1); i <= max_x; i = i + 1 ) {
        
        indices = find(x == i);
        
        yi = y.elem(indices) ;
        pi = p.elem(indices) ;
        
        for (int j = 0; j < yi.n_elem; j = j + 1 ) {
            for (int k = 0; k < non_zero_ys.n_elem; k = k + 1 ) {
                if (non_zero_ys[k] + yi[j] <= threshold) {
                    
                    // take existing probability for current y, and add new value
                    p_tmp[non_zero_ys[k] + yi[j]] = p_tmp[non_zero_ys[k] + yi[j]] + (p_out[non_zero_ys[k]] * pi[j]) ;
                    
                } else {
                    
                    // if value greater than threshold, then place in last slot of output vector
                    p_tmp[threshold] = p_tmp[threshold] + (p_out[non_zero_ys[k]] * pi[j]) ;
                    
                }
            }
            
        }
        
        // update output vector p_out
        p_out = p_tmp ;
        // refill p_tmp with zeros
        p_tmp.zeros() ;
        // find the non-zero entries of p_out and the corresponding integer values (scores)
        non_zero_ps = find(p_out != 0) ;
        non_zero_ys = integers(non_zero_ps) ;
        
    }
    
    
    // only keep non-zero entries
    p_out = p_out(non_zero_ps) ;
    
    // print to console
    //cout << "the p vector is : " << p_out << endl ;
    
    // results are placed in matrix for nice output
    arma::mat output( p_out.n_elem, 0 ) ;
    output.insert_cols( 0, non_zero_ys ) ;
    output.insert_cols( 1, p_out ) ;
    
    // output is returned
    return output ;
}



