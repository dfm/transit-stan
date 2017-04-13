#ifndef TRANSIT_STAN_LIMB_DARKENING_HPP
#define TRANSIT_STAN_LIMB_DARKENING_HPP

#define PAL_AB  T a = (p-z)*(p-z), \
                       b = (p+z)*(p+z);
#define PAL_K   T k0=acos((p2+z2-1.0)/(2.0*p*z)), \
                  k1=acos((1.0+z2-p2)/(2.0*z));

#define PAL_CI  T ci = 2.0/(9.0*M_PI*sqrt(1.0-a)), \
                  cik = (1.0-5.0*z2+p2+a*b), \
                  cie = (z2+7.0*p2-4.0)*(1.0-a), \
                  cip = -3.0*(p+z)/(p-z); \
                fk = ci*cik; fe = ci*cie; fp = ci*cip; \
                kk = ee = pp = 1; \
                k = sqrt(4.0*z*p/(1.0-a)); \
                n = (a-b)/a; \
                f0 = p2; \
                f2 = 0.5*f0*(f0+2.0*z2);

#define PAL_CG  PAL_K \
                T cg = 1.0/(9.0*M_PI*sqrt(p*z)), \
                       cgk = ((1.0-b)*(2.0*b+a-3.0)-3.0*(p+z)*(p-z)*(b-2.0)), \
                       cge = 4.0*p*z*(z2+7.0*p2-4.0), \
                       cgp = -3.0*(p+z)/(p-z); \
                fk = cg*cgk; fe = cg*cge; fp = cg*cgp; \
                kk = ee = pp = 1; \
                k = sqrt((1.0-a)/(4.0*p*z)); \
                n = (a-1.0)/a; \
                f0 = (p2*k0+k1-sqrt(z2-0.25*(1.0+z2-p2)*(1.0+z2-p2)))/M_PI; \
                f2 = (k1+p2*(p2+2.0*z2)*k0-0.25*(1.0+5.0*p2+z2)*sqrt((1.0-a)*(b-1.0)))/(2.0*M_PI);

template <typename T1, typename T2, typename T3, typename T4>
typename boost::math::tools::promote_args<T1, T2, T3, T4>::type
limb_darkening(const T1& q1_0, const T2& q2_0, const T3& p, const T4& z0, std::ostream* pstream__) {
  using std::abs;
  using std::sqrt;
  using std::acos;
  typedef typename boost::math::tools::promote_args<T1, T2, T3, T4>::type T;

  int kk = 0, ee = 0, pp = 0;
  T1 q1 = sqrt(q1_0);
  T2 q2 = 2.0 * q2_0;
  T u1 = q1 * q2, u2 = q1 * (1.0 - q2),
    w0, w1, w2,
    f0 = T(0.0), f1 = T(0.0), fk = T(0.0), fe = T(0.0), fp = T(0.0),
    f2 = T(0.0), n = T(0.0), k = T(0.0),
    df, z2, p2;

  // Pre-compute the constant coefficients.
  w0 = (6.0-6.0*u1-12.0*u2) / (6.0-2.0*u1-u2);
  w1 = ( 6.0*u1+12.0*u2   ) / (6.0-2.0*u1-u2);
  w2 = (   6.0*u2         ) / (6.0-2.0*u1-u2);
  T z = abs(z0);
  z2 = z * z;
  p2 = p * p;

  // Run through the cases and compute the coefficients from Table A1&2.
  if (z <= 0.0 && p <= 1.0) {
    // M&A 10, Pal A
    f0 = p2;
    f1 = 2.0/3.0*(1.0-(1.0-f0)*sqrt(1.0-f0));
    f2 = 0.5*f0*f0;
  } else if (z <= p-1.0) {
    // M&A 11, Pal A_G
    f0 = T(1.0);
    f1 = T(2.0/3.0);
    f2 = T(0.5);
  } else if (z < p && z < 1.0 - p - DBL_EPSILON) {
    // M&A 9, Pal B
    PAL_AB PAL_CI
      f1 = T(2.0/3.0);
  } else if (z < p && abs(z-1.0+p) <= DBL_EPSILON) {
    // M&A -, Pal B_T
    f0 = p2;
    f1 = (2.0/(3.0*M_PI)*acos(1.0-2.0*p)-4.0/(9.0*M_PI)
        * (3.0+2.0*p-8.0*f0)*sqrt(p*(1.0-p)));
    f2 = 0.5*f0*(f0+2.0*z2);
  } else if (z < p) {
    // M&A 8, Pal B_G
    PAL_AB PAL_CG
      f1 = T(2.0/3.0);
  } else if (abs(z-p) <= DBL_EPSILON && z < 1.0-p-DBL_EPSILON) {
    // M&A 5, Pal C
    T t = T(2.0/(9.0*M_PI));
    f0 = p2;
    f1 = T(1.0/3.0);
    fk = t*(1.0-4.0*f0); kk = 1;
    fe = t*4.0*(2.0*f0-1.0); ee = 1;
    f2 = 1.5*f0*f0;
    k = 2.0*p;
  } else if (abs(z-p) <= DBL_EPSILON && abs(z-1.0+p) <= DBL_EPSILON) {
    // M&A 6, Pal C_T
    f0 = T(0.25);
    f1 = T(1.0/3.0-4.0/(9.0*M_PI));
    f2 = T(3.0/32.0);
  } else if (abs(z-p) <= DBL_EPSILON) {
    // M&A 7, Pal C_G
    PAL_AB PAL_K
      f0 = (p2*k0+k1-sqrt(z2-0.25*(1.0+z2-p2)*(1.0+z2-p2)))/M_PI;
    f1 = T(1.0/3.0);
    fk = -(1.0-4.0*p2)*(3.0-8.0*p2)/(9.0*M_PI*p); kk = 1;
    fe = 16.0*p*(2.0*p2-1.0)/(9.0*M_PI); ee = 1;
    f2 = (k1+p2*(p2+2.0*z2)*k0-0.25*(1.0+5.0*p2+z2)*sqrt((1.0-a)*(b-1.0)))/(2.0*M_PI);
    k = 1.0/(2.0*p);
  } else if (z < 1.0-p-DBL_EPSILON) {
    // M&A 3, Pal D
    PAL_AB PAL_CI
  } else if (abs(z-1.0+p) <= DBL_EPSILON) {
    // M&A 4, Pal E
    f0 = p2;
    f1 = (2.0/(3.0*M_PI)*acos(1.0-2.0*p)-4.0/(9.0*M_PI)
        *(3.0+2.0*p-8.0*f0)*sqrt(p*(1.0-p)));
    f2 = 0.5*f0*(f0+2.0*z2);
  } else if (z < 1.0+p-DBL_EPSILON) {
    // M&A 2, Pal F
    PAL_AB PAL_CG
  }

  // Compute the base flux change.
  df = w0*f0+w1*f1+w2*f2;

  // Add in the elliptic integral terms.
  if (abs(w1) > DBL_EPSILON) {
    if (kk && abs(fk) > DBL_EPSILON)
      df += w1 * fk * ellint_1(k, pstream__);
    if (ee && abs(fe) > DBL_EPSILON)
      df += w1 * fe * ellint_2(k, pstream__);
    if (pp && abs(fp) > DBL_EPSILON)
      df += w1 * fp * ellint_3(-k, n, pstream__);
  }

  return 1.0 - df;
}

#undef PAL_AB
#undef PAL_K
#undef PAL_CI
#undef PAL_CG

#endif
