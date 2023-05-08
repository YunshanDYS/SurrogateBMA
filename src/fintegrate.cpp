// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppNumerical.h> 

using namespace Numer;

using Rcpp::NumericVector;

class Mintegrand: public Func
{
private:
    const int m;
    const double mu;
    const double sig2;
public:
    Mintegrand(int m_, double mu_, double sig2_) : m(m_), mu(mu_), sig2(sig2_) {}

    double operator()(const double& x) const
    {
        double gx;
        double px;
        if (m == 4) {
            gx = std::sqrt(x);
        } else {
            // m == 5
            gx = std::log(x + 1);
        }
        px = R::dnorm(x, mu, std::sqrt(sig2), FALSE);
        const double f = gx * px;
        return f;
    }
};

// [[Rcpp::export]]
double gpxi_int(int i, int m, NumericVector beta, NumericVector S0, NumericVector S1, double sig2)
{
    int n0 = S0.length();
    int n1 = S1.length();
    double higher_tol = 50.0;
    double lower_tol;
    if (m == 4) {
        lower_tol = 0.0;
    } else if (m == 5) {
        lower_tol = -1.0;
    }
    double out;
    double mu;
    if (i == 1) {
        for (int j = 0; j < n1; j++) {
            mu = beta[0] + beta[1] + (beta[2] + beta[3]) * S1[j];
            Mintegrand f(m, mu, sig2);
            double err_est;
            int err_code;
            double fout = integrate(f, lower_tol, higher_tol, err_est, err_code);
            out = out + fout;
        }
        out = out / n1;
    } else if (i == 0) {
        for (int j = 0; j < n0; j++) {
            mu = beta[0] + beta[2] * S0[j];
            Mintegrand f(m, mu, sig2);
            double err_est;
            int err_code;
            double fout = integrate(f, lower_tol, higher_tol, err_est, err_code);
            out = out + fout;
        }
        out = out / n0;
    } else {
        // i == 2
        for (int j = 0; j < n0; j++) {
            mu = beta[0] + beta[1] + (beta[2] + beta[3]) * S0[j];
            Mintegrand f(m, mu, sig2);
            double err_est;
            int err_code;
            double fout = integrate(f, lower_tol, higher_tol, err_est, err_code);
            out = out + fout;
        }
        out = out / n0;
    }
    return out;
}
