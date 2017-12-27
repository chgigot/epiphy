#include <Rcpp.h>
using namespace Rcpp;

const double PI_ = 3.14159265358979323846;

// [[Rcpp::export]]
NumericVector ntotfoci(int nfoci, int xrate, int ngen, int ngen_active) {
    // Convention: if ngen_active == 0, all the foci are active at every generation.
    // Si ngen_active > ngen => error
    // Si ngen_active is negative => error
    NumericVector nfoci_new(ngen + 1);
    int nfoci_active(0);
    if (ngen_active == 0) {
        ngen_active = ngen;
    }
    nfoci_new[0] = nfoci;
    for (int i1 = 1; i1 <= ngen; ++i1) {
        for (int i2 = i1 - ngen_active; i2 <= i1; ++i2) {
            if (i2 >= 0) {
                nfoci_active += nfoci_new[i2];
            }
        }
        nfoci_new[i1] = nfoci_active * xrate;
        nfoci_active = 0;
    }
    return nfoci_new;
}

// [[Rcpp::export]]
NumericMatrix dispersalCPP(int nfoci, int xrate, double lambda, int ngen,
                           int ngen_active = 0) {
    if (ngen_active == 0) {
        ngen_active = ngen;
    }
    int nfoci_tot(0);
    NumericVector nfoci_new(ngen + 1);
    nfoci_new = ntotfoci(nfoci, xrate, ngen, ngen_active);
    for (int i1 = 0; i1 <= ngen; ++i1) {
        nfoci_tot += nfoci_new[i1];
    }
    NumericMatrix foci(nfoci_tot, 3);
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

    for (int i1 = 1; i1 <= ngen; ++i1) {
        NumericMatrix newfoci(xrate, 3);
        int ibeg(i1 - ngen_active);
        if (ibeg < 0) ibeg = 0;
        int beg(0);
        for (int i2 = 0; i2 < ibeg; ++i2) {
            beg += nfoci_new[i2];
        }
        int end(beg);
        for (int i2 = ibeg; i2 < i1; ++i2) {
            end += nfoci_new[i2];
        }
        //printf("[%d, %d[... beg = %d and end = %d\n", ibeg, i1, beg, end);
        for (int i2 = beg; i2 < end; ++i2) {
            focus   = foci(i2, _);
            newdir  = runif(xrate, 0, 2*PI_);
            newdist = rexp(xrate, lambda);
            newfoci(_, 0) = rep(focus[0], xrate) + newdist * cos(newdir);
            newfoci(_, 1) = rep(focus[1], xrate) + newdist * sin(newdir);
            newfoci(_, 2) = rep(i1, xrate);
            int i3_base(end + (i2 - beg) * xrate);
            //printf("--> base : %d, to: %d\n", i3_base, (end + (i2 + 1) * xrate));
            for (int i3 = i3_base;
                     i3 < (end + (i2 - beg + 1) * xrate);
                     i3++) {
                //printf("i3: %d\n", i3);
                //printf("i3: %d, nr: %d... i3-i3d: %d, nr: %d\n",
                //       i3, foci.nrow(), i3 - i3_base,  newfoci.nrow());
                //printf("0\n");
                foci(i3, 0) = newfoci(i3 - i3_base, 0);
                //printf("1\n");
                foci(i3, 1) = newfoci(i3 - i3_base, 1);
                //printf("2\n");
                foci(i3, 2) = newfoci(i3 - i3_base, 2);
                //printf("3\n");
            }
        }
    }
    return foci;
}

// [[Rcpp::export]]
NumericMatrix collectCPP(NumericMatrix disperse, NumericMatrix quadrat) {
    double count(0.0), incidence(0.0);
    LogicalVector inx(disperse.nrow()), iny(disperse.nrow());
    NumericMatrix res(quadrat.nrow(), 8);

    for (int i1 = 0; i1 < quadrat.nrow(); ++i1) {
        inx = (quadrat(i1, 2) < disperse(_, 0)) & (disperse(_, 0) <= quadrat(i1, 4));
        iny = (quadrat(i1, 3) < disperse(_, 1)) & (disperse(_, 1) <= quadrat(i1, 5));
        count = sum(inx * iny); // Nice trick!
        if (count == 0.0) {
            incidence = 0.0;
        } else {
            incidence = 1.0;
        }
        for (int i2 = 0; i2 < 6; ++i2) {
            res(i1, i2) = quadrat(i1, i2);
        }
        res(i1, 6) = count;
        res(i1, 7) = incidence;
    }
    return res;
}






















