#include <Rcpp.h>
#include <cmath>

// Square supported kernel:
double kern(const Rcpp::NumericVector & site, const double constant = 20.3/4.0) { // # 20.3 # not good (aire sous la courb !!!!)
    double exponent(0.0);
    if ((-1 <= site(0)) && (site(0) <= 1) && // x coord.
        (-1 <= site(1)) && (site(1) <= 1)) { // y coord.
        //exponent = (1/(1 - std::pow(site(0), 2)) +
        //    1/(1 - std::pow(site(1), 2)));
        exponent = (1/(1 - (site(0) * site(0))) +
            1/(1 - (site(1) * site(1))));
        return constant * std::exp(-exponent);
    } else {
        return 0.0;
    }
}

// h-scaled and renormalized kernel:
double kern_hscaled(const Rcpp::NumericVector & node,
                    const Rcpp::NumericVector & site,
                    const double bandwidth, const double denom) {
    Rcpp::NumericVector rdist(2);
    for (int i1 = 0; i1 < 2; ++i1) {
        rdist(i1) = (node(i1) - site(i1)) / bandwidth;
    }
    if (denom == 0) {
        // Without edge correction:
        return kern(rdist) / (bandwidth * bandwidth);
    } else {
        // With edge correction:
        return kern(rdist) / denom;
    }
}

// [[Rcpp::export]]
Rcpp::NumericVector p_hscaled(const Rcpp::DataFrame & nodes,
                              const Rcpp::DataFrame & sites,
                              double bandwidth, bool edgeCorrection = false) {

    int nRowNodes(nodes.nrow()), nRowSites(sites.nrow());
    Rcpp::NumericVector xNodes = nodes(0), yNodes = nodes(1);
    Rcpp::NumericVector xSites = sites(0), ySites = sites(1), iSites = sites(2);
    double subDenom(0.0), numer(0.0);
    Rcpp::NumericVector res(nRowNodes);
    // We only use the two first element of a node/site, i.e. the
    // coordinates x and y:
    Rcpp::NumericVector tmpNode(2), tmpSite(2), rdist(2);
    for (int i1 = 0; i1 < nRowNodes; ++i1) {
        numer = 0.0;
        for (int i2 = 0; i2 < nRowSites; ++i2) {
            tmpNode(0) = xNodes(i1);
            tmpNode(1) = yNodes(i1);
            tmpSite(0) = xSites(i2);
            tmpSite(1) = ySites(i2);
            //for (int i3 = 0; i3 < 2; ++i3) {
            //    tmpNode(i3) = nodes.at(i1, i3);
            //    tmpSite(i3) = sites.at(i2, i3);
            //}
            subDenom = 0.0;
            if (edgeCorrection) {
                for (int i4 = 0; i4 < nRowNodes; ++i4) {
                    rdist(0) = (xNodes(i4) - tmpSite(0)) / bandwidth;
                    rdist(1) = (yNodes(i4) - tmpSite(1)) / bandwidth;
                    //for (int i3 = 0; i3 < 2; ++i3) {
                    //    rdist(i3) = (nodes.at(i4, i3) - tmpSite(i3)) / bandwidth;
                    //}
                    subDenom += kern(rdist);
                }
            }
            numer += kern_hscaled(tmpNode, tmpSite, bandwidth, subDenom) * iSites(i2);// sites.at(i2, 2);
        }
        double denom = Rcpp::sum(iSites);//Rcpp::sum(sites(Rcpp::_, 2));
        res(i1) = numer / denom;
    }
    return res;
}




