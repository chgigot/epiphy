#include <Rcpp.h>
using namespace Rcpp;

const double PI_ = 3.14159265358979323846;

int ntotfoci(int nfoci, int xrate, int ngen) {
    while (ngen > 0) {
        nfoci = nfoci + nfoci * xrate;
        --ngen;
    }
    return nfoci;
}

// [[Rcpp::export]]
NumericMatrix dispersalCPP(int nfoci, int xrate, int ngen, double lambda) {

    long nrow(ntotfoci(nfoci, xrate, ngen));
    NumericMatrix foci(nrow, 3);
    NumericMatrix initfoci(nfoci, 3);
    NumericVector focus(3);
    NumericVector newdir(xrate), newdist(xrate);

    initfoci(_, 0) = runif(nfoci, 0, 1);
    initfoci(_, 1) = runif(nfoci, 0, 1);
    initfoci(_, 2) = rep(0, nfoci);
    for (int i1 = 0; i1 < nfoci; i1++) {
        foci(i1, 0) = initfoci(i1, 0);
        foci(i1, 1) = initfoci(i1, 1);
        foci(i1, 2) = initfoci(i1, 2);
    }

    for (int i1 = 0; i1 < ngen; i1++) {
        NumericMatrix newfoci(nfoci * xrate, 3);
        for (int i2 = 0; i2 < nfoci; i2++) {
            focus   = foci(i2, _);
            newdir  = runif(xrate, 0, 2*PI_);
            newdist = rexp(xrate, lambda);
            newfoci(_, 0) = rep(focus[0], xrate) + newdist * cos(newdir);
            newfoci(_, 1) = rep(focus[1], xrate) + newdist * sin(newdir);
            newfoci(_, 2) = rep(i1 + 1, xrate);
            int i3_base(nfoci + i2 * xrate);
            for (int i3 = i3_base;
                     i3 < (nfoci + (i2 + 1) * xrate);
                     i3++) {
                foci(i3, 0) = newfoci(i3 - i3_base, 0);
                foci(i3, 1) = newfoci(i3 - i3_base, 1);
                foci(i3, 2) = newfoci(i3 - i3_base, 2);
            }
        }
        nfoci = nfoci + nfoci * xrate;
    }
    return foci;
}

// NumericMatrix collectCPP(NumericMatrix disperse, NumericMatrix quadrat) {
//     NumericVector quad(6);
//     double count(0.0), incidence(0.0);
//     LogicalVector inx(disperse.nrow()), iny(disperse.nrow());
//     NumericMatrix res(quadrat.nrow(), 8);
//     for (int i1 = 0; i1 < quadrat.nrow(); ++i1) {
//         quad = quadrat(i1, _);
//         inx = (quad[2] < disperse(_, 0)) & (disperse(_, 0) <= quad[4]);
//         iny = (quad[3] < disperse(_, 1)) & (disperse(_, 1) <= quad[5]);
//         count = sum(inx * iny); // Nice trick!
//         if (count == 0.0) {
//             incidence = 0.0;
//         } else {
//             incidence = 1.0;
//         }
//         for (int i2 = 0; i2 < 6; ++i2) {
//             res(i1, i2) = quad[i2];
//         }
//         res(i1, 6) = count;
//         res(i1, 7) = incidence;
//     }
//     return res;
// }






















