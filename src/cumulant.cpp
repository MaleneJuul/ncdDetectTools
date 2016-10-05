#include <Rcpp.h>
using namespace Rcpp;

double c0(double t, double p, double s){
  return log(1-p+p*exp(s*t));
}

double c1(double t, double p, double s){
  double num = s*p*exp(s*t);
  double den = 1-p+p*exp(s*t);
  return num/den;
}

double c2(double t, double p, double s){
  double num = s*s*p*exp(s*t)*(1-p);
  double sqrtden = (1-p+p*exp(s*t));
  return num/sqrtden/sqrtden;
}

// [[Rcpp::export]]
double cumulant(double t, NumericVector p, NumericVector s) {
  /* Return cumulant generating function*/
  double ret = 0;
  for(int i = 0; i < p.size(); ++i)
    ret += c0(t, p[i], s[i]);
  return ret;
}

// [[Rcpp::export]]
double cumulantD1(double t, NumericVector p, NumericVector s) {
  /* Return first derivative of cumulant generating function */
  double ret = 0;
  for(int i = 0; i < p.size(); ++i)
    ret += c1(t, p[i], s[i]);
  return ret;
}

// [[Rcpp::export]]
double cumulantD2(double t, NumericVector p, NumericVector s) {
  /* Return second derivative of cumulant generating function */
  double ret = 0;
  for(int i = 0; i < p.size(); ++i)
    ret += c2(t, p[i], s[i]);
  return ret;
}

