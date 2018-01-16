#include <Rcpp.h>
#include <cstdlib>
#include <string>
using namespace Rcpp;

// [[Rcpp::export]]
double costTotiCPP(long i1, NumericMatrix & flow, NumericMatrix & cost,
                   std::string type = "both", bool average = false) {
    if (!((type == "in") || (type == "out") || (type == "both"))) {
        exit(EXIT_FAILURE);
    }
    i1 -= 1; // to deal with the fact taht in C++, indices start at 0 (and not in R)
    double res(0);
    //    sorties             entrées
    if ((sum(flow(i1, _)) <= sum(flow(_, i1))) && ((type == "in") || (type == "both"))) {
        res = -sum(flow(_, i1) * cost(_, i1));
        if (average) res /= sum(flow(_, i1)); // sûr que i1 est à droite?
        //       sorties            entrées
    } else if ((sum(flow(i1, _)) > sum(flow(_, i1))) && ((type == "out") || (type == "both"))) {
        res = sum(flow(i1, _) * cost(i1, _));
        if (average) res /= sum(flow(i1, _));
    }
    return res; // "can return 0 if no if verified"
}

// [[Rcpp::export]]
double costTotCPP(NumericMatrix & flow, NumericMatrix & cost) {
    long res(0);
    for (long i1 = 1; i1 <= flow.nrow(); ++i1) {
        res += costTotiCPP(i1, flow, cost, "in");
    }
    return res;
}

// [[Rcpp::export]]
NumericMatrix optTransDT2MatCPP(NumericMatrix x, int suLen) {
    NumericMatrix res(suLen, suLen); // By default, all elements set to 0.0 (http://www.thecoatlessprofessor.com/programming/rcpp/unofficial-rcpp-api-docs/)
    for (int i1 = 0; i1 < x.nrow(); ++i1) {
        std::cout << i1 << std::endl;
        //std::cout << x(i1, 0) << " " << x(i1, 1) << " = " << x(i1, 2) << std::endl;
        res(x(i1, 0), x(i1, 1)) = (double) x(i1, 2);
    }
    return res;
}














