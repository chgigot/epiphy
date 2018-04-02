#include <Rcpp.h>
#include <string>

// [[Rcpp::export]]
double costTotiCPP(long i1, const Rcpp::NumericMatrix & flow,
                   const Rcpp::NumericMatrix & cost,
                   bool averaged = true,    // TODO: used?
                   bool absolute = false) { // TODO: used?
    i1 -= 1; // Because indices start at 0 in C++, and 1 in R.
    double res(0), inflow(0.0), outflow(0.0), sign(1.0), denom(1.0);
    // A donor unit (outflow > 0):
    if ((outflow = Rcpp::sum(flow(i1, Rcpp::_))) > 0) {
        denom = (averaged ? outflow : 1.0);
        res += Rcpp::sum(flow(i1, Rcpp::_) * cost(i1, Rcpp::_)) / denom;
    }
    // A receiver unit (inflow > 0):
    if ((inflow = Rcpp::sum(flow(Rcpp::_, i1))) > 0) {
        sign  = (absolute ? 1.0 : -1.0);
        denom = (averaged ? inflow : 1.0);
        res += sign * Rcpp::sum(flow(Rcpp::_, i1) * cost(Rcpp::_, i1)) / denom;
    }
    return res;
}

// [[Rcpp::export]]
Rcpp::NumericVector costTotCPP(const Rcpp::NumericMatrix & flow,
                               const Rcpp::NumericMatrix & cost) {
    Rcpp::NumericVector res(flow.nrow());
    for (int i1 = 0; i1 < flow.nrow(); ++i1) {
        res.at(i1) = costTotiCPP(i1 + 1, flow, cost);
        // +1 because it is here a R-like index.
    }
    return res;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix as_matrix_transport(Rcpp::List & x, double dim_mat) {
    Rcpp::NumericMatrix res(dim_mat); // A square matrix.
    // By default, a matrix constructed from dimensions will always be
    // initialized with all entries being zero. So no need to use res.fill(0.0);
    Rcpp::IntegerVector from(x[0]), to(x[1]);
    Rcpp::NumericVector mass(x[2]);
    for (int i1 = 0; i1 < from.length(); ++i1) {
        res.at((from.at(i1) - 1), (to.at(i1) - 1)) = mass.at(i1);
    }
    for (int i1 = 0; i1 < dim_mat; ++i1) {
        res(i1, i1) = 0; // No flows on the diagonal, i.e. not flows between a
        // unit and itself.
    }
    return res;
}
















