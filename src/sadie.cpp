#include <Rcpp.h>
#include <cstdlib>
#include <string>
using namespace Rcpp;

using namespace std;
// Utiliser costToti(...., "both") partout, pour prendre tous les flux en compte !!!

// [[Rcpp::export]]
double costTotiCPP(long i1, NumericMatrix & flow, NumericMatrix & cost,
                   std::string type = "both", bool averaged = true,
                   bool absolute = false) { // TODO: type not used any more.
    i1 -= 1; // to deal with the fact taht in C++, indices start at 0 (and not in R)
    flow(i1, i1) = 0; // Nothing for the cell itself (it not a flow)... pour la division, sinon erreur
    double res(0);
    // A receiver unit:
    // Inflow:
    if (sum(flow(_, i1)) > 0) {
        if (absolute) {
            if (averaged) {
                res += sum(flow(_, i1) * cost(_, i1)) / sum(flow(_, i1));
            } else {
                res += sum(flow(_, i1) * cost(_, i1));
            }
        } else {
            if (averaged) {
                res += -sum(flow(_, i1) * cost(_, i1)) / sum(flow(_, i1));
            } else {
                res += -sum(flow(_, i1) * cost(_, i1));
            }
        }
    }
    // A donor unit:
    // Outflow:
    if (sum(flow(i1, _)) > 0) {
        if (averaged) {
            res += sum(flow(i1, _) * cost(i1, _)) / sum(flow(i1, _));
        } else {
            res += sum(flow(i1, _) * cost(i1, _));
        }
    }
    return res; // "can return 0 if no if verified"
}

// [[Rcpp::export]]
double costTotCPP(NumericMatrix & flow, NumericMatrix & cost) {
    long res(0);
    for (long i1 = 1; i1 <= flow.nrow(); ++i1) {
        res += costTotiCPP(i1, flow, cost, "donor");
    }
    return res;
}

// [[Rcpp::export]]
NumericMatrix as_matrix_transport(List & x, double dim_mat) {
    NumericMatrix res(dim_mat); // By default, the matrices constructed from dimensions will always be initialized with all entries being zero (0)
    //res.fill(0.0); // Create a square matrix (all elements set to 0.0)
    IntegerVector from(x[0]), to(x[1]);
    NumericVector mass(x[2]);
    for (int i1 = 0; i1 < from.length(); ++i1) {
        res.at((from.at(i1) - 1), (to.at(i1) - 1)) = mass.at(i1);
    }
    return res;
}
















