#ifndef TRANSIT_STAN_ELLINT_HPP
#define TRANSIT_STAN_ELLINT_HPP

#include <cmath>
#include <stan/math/rev/scal.hpp>

#define ELLINT_CONV_TOL 1.0e-8
#define ELLINT_MAX_ITER 200

// K: 1.0 - k^2 >= 0.0
inline double ellint_1 (double k, std::ostream* pstream__) {
  if (1.0 < k*k) domain_error("ellint_1", "k is", k, "");
  double kc = std::sqrt(1.0 - k * k), m = 1.0, h;
  for (int i = 0; i < ELLINT_MAX_ITER; ++i) {
    h = m;
    m += kc;
    if (std::abs(h - kc) / h <= ELLINT_CONV_TOL) break;
    kc = std::sqrt(h * kc);
    m *= 0.5;
  }
  //if (i >= ELLINT_MAX_ITER)
  //  std::cerr << "Warning: ellint_1 did not converge; use result with caution\n";
  return M_PI / m;
}

// E: 1.0 - k^2 >= 0.0
inline double ellint_2 (const double k, std::ostream* pstream__) {
  if (1.0 < k*k) domain_error("ellint_2", "k is", k, "");
  double b = 1.0 - k * k, kc = std::sqrt(b), m = 1.0, c = 1.0, a = b + 1.0, m0;
  for (int i = 0; i < ELLINT_MAX_ITER; ++i) {
    b = 2.0 * (c * kc + b);
    c = a;
    m0 = m;
    m += kc;
    a += b / m;
    if (std::abs(m0 - kc) / m0 <= ELLINT_CONV_TOL) break;
    kc = 2.0 * std::sqrt(kc * m0);
  }
  //if (i >= ELLINT_MAX_ITER)
  //  std::cerr << "Warning: ellint_2 did not converge; use result with caution\n";
  return M_PI_4 * a / m;
}

// Pi: 1.0 - k^2 >= 0.0 & n < 1.0
inline double ellint_3 (double k, double n, std::ostream* pstream__) {
  if (1.0 < k*k) domain_error("ellint_3", "k is", k, "");
  if (n >= 1.0) domain_error("ellint_3", "n is", n, "");
  double kc = std::sqrt(1.0 - k * k), p = std::sqrt(1.0 - n), m0 = 1.0, c = 1.0,
          d = 1.0 / p, e = kc, f, g;
  for (int i = 0; i < ELLINT_MAX_ITER; ++i) {
    f = c;
    c += d / p;
    g = e / p;
    d = 2.0 * (f * g + d);
    p = g + p;
    g = m0;
    m0 = kc + m0;
    if (std::abs(1.0 - kc / g) <= ELLINT_CONV_TOL) break;
    kc = 2.0 * std::sqrt(e);
    e = kc * m0;
  }
  //if (i >= ELLINT_MAX_ITER)
  //  std::cerr << "Warning: ellint_3 did not converge; use result with caution\n";
  return M_PI_2 * (c * m0 + d) / (m0 * (m0 + p));
}

#undef ELLINT_CONV_TOL
#undef ELLINT_MAX_ITER

inline var ellint_1 (const var& z_var, std::ostream* pstream__) {
  double z = z_var.val();
  double Kz = ellint_1(z, pstream__), Ez = ellint_2(z, pstream__);
  return var(new precomp_v_vari(
    Kz, z_var.vi_, (Ez / (1.0 - z*z) - Kz) / z
  ));
}

inline var ellint_2 (const var& z_var, std::ostream* pstream__) {
  double z = z_var.val();
  double Kz = ellint_1(z, pstream__), Ez = ellint_2(z, pstream__);
  return var(new precomp_v_vari(
    Ez, z_var.vi_, (Ez - Kz) / z
  ));
}

inline var ellint_3 (const var& k_var, const var& n_var, std::ostream* pstream__) {
  double k = k_var.val();
  double n = n_var.val();
  double Kk = ellint_1(k, pstream__), Ek = ellint_2(k, pstream__), Pnk = ellint_3(k, n, pstream__),
          k2 = k * k, n2 = n * n;
  return var(new precomp_vv_vari(
    Pnk, k_var.vi_, n_var.vi_,
    -k * (Ek / (k2 - 1.0) + Pnk) / (k2-n),
    (0.5*(Ek + (Kk*(k2-n) + Pnk*(n2-k2))/n)/(n-1.0))/(k2-n)
  ));
}

inline var ellint_3 (double k, const var& n_var, std::ostream* pstream__) {
  double n = n_var.val();
  double Kk = ellint_1(k, pstream__), Ek = ellint_2(k, pstream__), Pnk = ellint_3(k, n, pstream__),
          k2 = k * k, n2 = n * n;
  return var(new precomp_v_vari(
    Pnk, n_var.vi_,
    (0.5*(Ek + (Kk*(k2-n) + Pnk*(n2-k2))/n)/(n-1.0))/(k2-n)
  ));
}

inline var ellint_3 (const var& k_var, double n, std::ostream* pstream__) {
  double k = k_var.val();
  double Ek = ellint_2(k, pstream__), Pnk = ellint_3(k, n, pstream__),
          k2 = k * k;
  return var(new precomp_v_vari(
    Pnk, k_var.vi_, -k * (Ek / (k2 - 1.0) + Pnk) / (k2-n)
  ));
}

#endif
