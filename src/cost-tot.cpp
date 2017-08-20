#include <Rcpp.h>
#include <cstdlib>
#include <string>
using namespace Rcpp;

// [[Rcpp::export]]
double costTotiCPP(long i, NumericMatrix & flow, NumericMatrix & cost,
                   std::string type = "both", bool average = false) {
    if (!((type == "in") || (type == "out") || (type == "both"))) {
        exit(EXIT_FAILURE);
    }
    i -= 1; // to deal with the fact taht in C++, indices start at 0 (and not in R)
    double res(0);
    //    sorties             entrées
    if ((sum(flow(i, _)) <= sum(flow(_, i))) && ((type == "in") || (type == "both"))) {
        res = -sum(flow(_, i) * cost(_, i));
        if (average) res /= sum(flow(_, i)); // sûr que i est à droite?
        //       sorties            entrées
    } else if ((sum(flow(i, _)) > sum(flow(_, i))) && ((type == "out") || (type == "both"))) {
        res = sum(flow(i, _) * cost(i, _));
        if (average) res /= sum(flow(i, _));
    }
    return res; // "can return 0 if no if verified"
}

// [[Rcpp::export]]
double costTotCPP(NumericMatrix & flow, NumericMatrix & cost) {
    long res(0);
    for (long i = 1; i <= flow.nrow(); ++i) {
        res += costTotiCPP(i, flow, cost, "in");
    }
    return res;
}

// [[Rcpp::export]]
NumericVector insertCPP(NumericVector & vect, int i, // i = position
                        const double & value) {
    i -= 1; // to deal with the fact taht in C++, indices start at 0 (and not in R)
    vect.insert(i, value); // ref: http://statr.me/rcpp-note/api/Vector_funs.html
    return vect;
}

