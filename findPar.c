// Compile:
//   gcc -O3 -march=native -flto -fopenmp -fno-math-errno -std=c11 findPar.c -lm -o rsd
//
// Run:
//   OMP_NUM_THREADS=16 ./rsd n k h [binFlag] [xiCoef]
//
// Example:
//   ./rsd 1024 652 106
//   OMP_NUM_THREADS=16 ./rsd 1024 652 106
//   OMP_NUM_THREADS=16 ./rsd 1024 652 106 1
//   OMP_NUM_THREADS=16 ./rsd 1024 652 106 2.5

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
    double best_log2_cost;

    int n_adj, k_adj, h, beta, numTrun;
    int omega, flin;

    int f1, mu1, eta1;
    int f2, mu2, eta2;
    int f_pre, nu, eta_pre;

    double p1, p2;

    int n1, n2;
    int64_t gamma;
    int64_t D1, d1, N1, m1;                 // Di and di stnad for Delta_{i} and delta_{i}.
    int64_t D2, d2, N2, N2p, N2pp, m2, m2p;
    
    double T11, T12, T2, T31, T32;
    double com1, com2, com3;
} BestResultSearchParams;

static inline int64_t choose2_i64(int64_t x) {
    if (x < 2) return 0;
    return (x * (x - 1)) / 2;
}

static inline void compute_w_flin(int n, int k, int h, int *w, int *flin) {
    *w = (n-k) / h;
    *flin = (n-k) % h;
}

static inline void compute_mu_f(int eta, int h, int *mu, int *f) {
    *mu = eta / h;
    *f  = eta % h;
}

static inline double p_f_mu(int f, int mu, int beta, int h) {
    double b = (double)beta;
    double a1 = 1.0 - ((double)(mu + 1)) / b;
    double a0 = 1.0 - ((double)mu) / b;
    if (a1 <= 0.0 || a0 <= 0.0) return 0.0;
    return pow(a1, (double)f) * pow(a0, (double)(h - f));
}

static inline int64_t compute_Delta(int h, int beta, int mu, int omega, int flin, int f) {
    int64_t t = (int64_t)beta - (int64_t)mu - (int64_t)omega - 1;
    if (f + flin <= h) {
        return (int64_t)h * choose2_i64(t) + (int64_t)(h - f - flin) * t;
    } else {
        return (int64_t)h * choose2_i64(t) + (int64_t)(h - f - flin) * (t - 1);
    }
}

static inline int64_t compute_m(int h, int beta, int f, int mu) {
    int64_t a = (int64_t)beta - (int64_t)mu - 1;
    int64_t b = (int64_t)beta - (int64_t)mu;
    return (int64_t)f * choose2_i64(a) + (int64_t)(h - f) * choose2_i64(b);
}

static inline double comGE(int64_t numRow, int64_t numCol, double constantGauss, int binFlag) {
    double t1, t2, t3;
    int64_t minVal = numRow;
    int64_t maxVal = numCol;
    if (minVal > numCol) {
        minVal = numCol;
        maxVal = numRow;
    }

    if (binFlag) {
        t1 = (double)numRow * (double)numCol * (double)minVal;
        t2 = pow(maxVal, constantGauss);
        t3 = (double)numRow * (double)numRow * (double)numCol / log2(numRow);
    }
    else {
        t1 = (double)numRow * (double)numCol * (double)minVal;
        t2 = pow(maxVal, constantGauss);
        t3 = t2; // t3 is not used in the non-binary case
    }

    if (t1 < t2) {
        if (t1 < t3) {
            return t1;
        }
        else {
            return t3;
        }
    }
    else {
        if (t2 < t3) {
            return t2;
        }
        else {
            return t3;
        }
    }
}

static inline int64_t xi(int64_t n, double xiVal, int mode) {
    // xi(N) function
    // mode == 2 for binary field     (it returns as xi(n) = xiVal * n)
    // mode == 3 for non-binary field (it returns as xi(n) = xiVal)

    if (mode == 1) {
        // [WWY+25]'s restriction
        return (int64_t)(1.5 * (double)n);
    }
    else if (mode == 2) {
        // For GF(2)
        return (int64_t)(xiVal * (double)n);
    }
    else if (mode == 3) {
        // For GF(2^{128})
        return xiVal;
    }
    else {
        // not defined mode.
        return 0;
    }
}

static inline void noGuessComplexity(int n, int k, int h, int beta, int numTrun, int omega, int flin, double constantGauss, BestResultSearchParams *best_global, double *best_cost, bool *found_any, int binFlag)
{
    int n0 = k;
    if (binFlag) n0 -= h;

    int64_t N0, m0, D0;

    N0 = choose2_i64((int64_t)n0) + 2LL * (int64_t)n0;
    if (binFlag) N0 -= (int64_t)n0;
    m0 = (int64_t)h * choose2_i64((int64_t)beta);
    D0 = compute_Delta(h, beta, 0, omega, flin, 0);

    if (m0 >= N0) {
        *found_any = true;
        double T1, T2;

        int HCol = n;
        if (binFlag) HCol -= h;
        T1 = comGE(n-k, HCol, constantGauss, binFlag);
        if (binFlag) T1 += (n-k) * n;
        T2 = comGE(N0-D0, N0-D0, constantGauss, binFlag); // it can find a solution with only submatrix (N0 x N0).

        *best_cost = log2(T1 + T2);
        best_global->best_log2_cost = *best_cost;

        best_global->n_adj = n;
        best_global->k_adj = k;
        best_global->h = h;
        best_global->beta = beta;
        best_global->numTrun = numTrun;
        best_global->omega = omega;
        best_global->flin = flin;

        best_global->f1 = 0;
        best_global->mu1 = 0;
        best_global->eta1 = 0;

        best_global->f2 = 0;
        best_global->mu2 = 0;
        best_global->eta2 = 0;

        best_global->f_pre = 0;
        best_global->nu = 0;
        best_global->eta_pre = 0;

        best_global->p1 = 1;
        best_global->p2 = 1;

        best_global->n1 = 0;
        best_global->n2 = 0;
        best_global->gamma = 0;

        best_global->D1 = 0;
        best_global->d1 = 0;
        best_global->N1 = 0;
        best_global->m1 = 0;

        best_global->D2 = 0;
        best_global->d2 = 0;
        best_global->N2 = 0;
        best_global->N2p = 0;
        best_global->N2pp = 0;
        best_global->m2 = 0;
        best_global->m2p = 0;

        best_global->T11 = 0;
        best_global->T12 = 0;
        best_global->T2  = 0;
        best_global->T31 = 0;
        best_global->T32 = 0;

        best_global->com1 = log2(T1);
        best_global->com2 = log2(T2);
        best_global->com3 = 0;
    }
}

int search_params_best_omp(
    int n_in,
    int k_in,
    int h,
    double constantGauss,
    int binFlag,
    double xiCoef,
    BestResultSearchParams *out_best
) {
    if (!out_best) return 0;
    if (h <= 0) return 0;

    int n = n_in;
    int k = k_in;
    int numTrun = 0;

    if (n % h != 0) {
        int beta0 = n / h; // floor
        numTrun = n - beta0 * h;
        n -= numTrun;
        k -= numTrun;
    }
    if (n <= 0 || k <= 0) return 0;

    int beta = n / h;
    if (beta < 2) return 0;

    int omega, flin;
    compute_w_flin(n, k, h, &omega, &flin);
    if (binFlag) omega += 1;

    // eta별 반복값 미리 계산(을 위한 공간 할당)
    int *mu_of_eta = (int *)malloc((size_t)(k + 1) * sizeof(int));
    int *f_of_eta  = (int *)malloc((size_t)(k + 1) * sizeof(int));
    double *p_of_eta = (double *)malloc((size_t)(k + 1) * sizeof(double));
    int64_t *D_of_eta = (int64_t *)malloc((size_t)(k + 1) * sizeof(int64_t));
    int64_t *m_of_eta = (int64_t *)malloc((size_t)(k + 1) * sizeof(int64_t));
    int64_t *gamma_pre = (int64_t *)malloc((size_t)(k + 1) * sizeof(int64_t));
    int64_t *removed_pre = (int64_t *)malloc((size_t)(k + 1) * sizeof(int64_t));

    if (!mu_of_eta || !f_of_eta || !p_of_eta || !D_of_eta || !m_of_eta || !gamma_pre || !removed_pre) {
        free(mu_of_eta); free(f_of_eta); free(p_of_eta); free(D_of_eta);
        free(m_of_eta); free(gamma_pre); free(removed_pre);
        return 0;
    }

    // eta별 반복값 미리 계산
    for (int eta = 0; eta <= k; ++eta) {
        int mu, f;
        compute_mu_f(eta, h, &mu, &f);
        mu_of_eta[eta] = mu;
        f_of_eta[eta] = f;
        p_of_eta[eta] = p_f_mu(f, mu, beta, h);
        D_of_eta[eta] = compute_Delta(h, beta, mu, omega, flin, f);
        m_of_eta[eta] = compute_m(h, beta, f, mu);

        // 아래는 eta_pre가 들어가야 되지만, 경우에 따른 값을 모두 계산한 테이블을 만드는 것
        if (binFlag) gamma_pre[eta] = choose2_i64((int64_t)eta);
        else gamma_pre[eta] = choose2_i64((int64_t)eta)+ (int64_t)eta;
        removed_pre[eta] = (int64_t)f * choose2_i64((int64_t)mu + 1) + (int64_t)(h - f) * choose2_i64((int64_t)mu);
    }

    BestResultSearchParams best_global;
    double best_cost = DBL_MAX;
    bool found_any = false;

    // 전수조사 안하는 것으로 가능 여부 확인 및 가능시 복잡도로 갱신 (안풀리면 갱신 x)
    noGuessComplexity(n,k,h,beta,numTrun,omega,flin,constantGauss,&best_global,&best_cost,&found_any,binFlag);

    // main job
    #pragma omp parallel
    {
        double best_local = DBL_MAX;
        BestResultSearchParams bestL;
        bool found_local = false;

        #pragma omp for collapse(2) schedule(dynamic)
        for (int eta_pre = 1; eta_pre <= k; ++eta_pre) {
            for (int eta = 0; eta <= k - 1; ++eta) {
                int nu = mu_of_eta[eta_pre];
                int f_pre = f_of_eta[eta_pre];
                int64_t gamma = gamma_pre[eta_pre];
                int64_t zeta = removed_pre[eta_pre];

                int mu2 = mu_of_eta[eta];
                int f2  = f_of_eta[eta];
                
                int blockMaxLen = (n-k) + eta + eta_pre;
                if (binFlag) blockMaxLen += h;
                blockMaxLen = (blockMaxLen + h-1) / h;  // ceiling 
                if (blockMaxLen > beta) continue;

                int64_t D2 = D_of_eta[eta];
                int64_t d2 = D2 - zeta;
                if (d2 < 0) continue;

                int n2 = k - eta;
                if (binFlag) n2 -= h;
                if (n2 <= 0 || n2 <= eta_pre) continue;
                if (n2 != eta_pre + 1) continue; // accelator for limited search

                int64_t N2 = choose2_i64((int64_t)n2) + 2LL * (int64_t)n2;
                if (binFlag) N2 -= (int64_t)n2;
                int64_t N2p = N2 - gamma - d2;
                int64_t N2pp = N2p - eta_pre;
                if (N2pp <= 0) continue;

                int64_t m2 = m_of_eta[eta];
                int64_t m2p = m2 - gamma - d2;
                int64_t m2pp = m2p - eta_pre;
                int64_t xiVal;
                //xiVal = xi(N2pp, 0, 1);
                if (binFlag) {
                    xiVal = xi(N2pp, xiCoef, 2);
                }
                else {
                    xiVal = xi(N2pp, xiCoef, 3);
                }
                int64_t hm2pp = N2pp + xiVal;
                if (m2pp <= 0) continue;

                if (n2 != eta_pre && m2pp - N2pp < xiVal) continue;

                double p2 = p_of_eta[eta];
                if (p2 <= 0.0) continue;

                double T32;
                T32 = comGE(hm2pp, N2pp, constantGauss, binFlag);
                
                if (!(T32 > 0.0) || !isfinite(T32)) continue;

                for (int eta1 = 0; eta1 <= eta; ++eta1) {
                    int mu1 = mu_of_eta[eta1];
                    int f1  = f_of_eta[eta1];
                    int eta2 = eta - eta1;
                    if (k <= eta1) continue;

                    double p1 = p_of_eta[eta1];
                    if (p1 <= 0.0) continue;

                    int64_t D1 = D_of_eta[eta1];
                    if (D1 < 0) continue;
                    int64_t d1 = D1 - zeta;

                    int n1 = k - eta1;
                    if (binFlag) n1 -= h;
                    if (n1 <= 0) continue;

                    int64_t N1 = choose2_i64((int64_t)n1) + 2LL * (int64_t)n1;
                    if (binFlag) N1 -= (int64_t)n1;
                    if (N1 <= 0) continue;

                    int64_t m1 = m_of_eta[eta1];
                    if (m1 <= 0) continue;

                    int64_t tau = (m1 - D1) - (m2 - D2);
                    if (tau < 0) continue;

                    int64_t hm2p = hm2pp + eta_pre;
                    int64_t hm1p = hm2p + tau;
                    int64_t hm1 = hm1p + gamma + d1;

                    double T11, T12, T2, T31, expected, cost_log2;

                    if (m1 >= N1) {
                        int HCol = n-eta1;
                        if (binFlag) HCol -= h;
                        T11 = comGE(n-k, HCol, constantGauss, binFlag);
                        if (binFlag) T11 += (n-k) * (n-eta1);
                        T12 = (double)D1 + (double)(m1 - D1) * (double)(n1+1) * (double)(n1+1);
                        T2 = comGE(m1 - D1, N1 - D1, constantGauss, binFlag);

                        expected = (T11 + T12 + T2) / p1;
                        cost_log2 = log2(expected);

                        if (cost_log2 < best_local) {
                            best_local = cost_log2;
                            found_local = true;

                            bestL.best_log2_cost = cost_log2;

                            bestL.n_adj = n;
                            bestL.k_adj = k;
                            bestL.h = h;
                            bestL.beta = beta;
                            bestL.numTrun = numTrun;
                            bestL.omega = omega;
                            bestL.flin = flin;

                            bestL.f1 = f1;
                            bestL.mu1 = mu1;
                            bestL.eta1 = eta1;

                            bestL.f2 = 0;
                            bestL.mu2 = 0;
                            bestL.eta2 = 0;

                            bestL.f_pre = 0;
                            bestL.nu = 0;
                            bestL.eta_pre = 0;

                            bestL.p1 = p1;
                            bestL.p2 = 1;

                            bestL.n1 = n1;
                            bestL.n2 = 0;
                            bestL.gamma = 0;

                            bestL.D1 = D1;
                            bestL.d1 = d1;
                            bestL.N1 = N1;
                            bestL.m1 = m1;

                            bestL.D2 = 0;
                            bestL.d2 = 0;
                            bestL.N2 = 0;
                            bestL.N2p = 0;
                            bestL.N2pp = 0;
                            bestL.m2 = 0;
                            bestL.m2p = 0;

                            bestL.T11 = T11;
                            bestL.T12 = T12;
                            bestL.T2  = T2;
                            bestL.T31 = 0;
                            bestL.T32 = 0;

                            bestL.com1 = log2((T11 + T12) / p1);
                            bestL.com2 = log2(T2 / p1);
                            bestL.com3 = 0;
                        }
                    }

                    int HCol = n-eta1;
                    if (binFlag) HCol -= h;
                    T11 = comGE(n-k, HCol, constantGauss, binFlag);
                    if (binFlag) T11 += (n-k) * (n-eta1);
                    T12 = (double)D1 + (double)(hm1 - D1) * (double)(n1+1) * (double)(n1+1);
                    T2 = comGE(hm1 - D1, N1 - D1, constantGauss, binFlag);
                    T31 = (double)hm2p * (double)(k - eta1 + 1) * (double)eta2;

                    double T1 = T11 + T12;
                    double T3 = T31 + T32;
                    expected = (T1 + T2) / p1 + T3 / p2;       
                    if (!(expected > 0.0) || !isfinite(expected)) continue;   

                    cost_log2 = log2(expected);
                    if (!isfinite(cost_log2)) continue;

                    if (cost_log2 < best_local) {
                        best_local = cost_log2;
                        found_local = true;

                        bestL.best_log2_cost = cost_log2;

                        bestL.n_adj = n;
                        bestL.k_adj = k;
                        bestL.h = h;
                        bestL.beta = beta;
                        bestL.numTrun = numTrun;
                        bestL.omega = omega;
                        bestL.flin = flin;

                        bestL.f1 = f1;
                        bestL.mu1 = mu1;
                        bestL.eta1 = eta1;

                        bestL.f2 = f2;      // f
                        bestL.mu2 = mu2;    // mu
                        bestL.eta2 = eta2;  // eta - eta1

                        bestL.f_pre = f_pre;
                        bestL.nu = nu;
                        bestL.eta_pre = eta_pre;

                        bestL.p1 = p1;
                        bestL.p2 = p2; // p total

                        bestL.n1 = n1;
                        bestL.n2 = n2;
                        bestL.gamma = gamma;

                        bestL.D1 = D1;
                        bestL.d1 = d1;
                        bestL.N1 = N1;
                        bestL.m1 = m1;

                        bestL.D2 = D2;
                        bestL.d2 = d2;
                        bestL.N2 = N2;
                        bestL.N2p = N2p;
                        bestL.N2pp = N2pp;
                        bestL.m2 = m2;
                        bestL.m2p = m2p;

                        bestL.T11 = T11;
                        bestL.T12 = T12;
                        bestL.T2  = T2;
                        bestL.T31 = T31;
                        bestL.T32 = T32;

                        bestL.com1 = log2(T1 / p1);
                        bestL.com2 = log2(T2 / p1);
                        bestL.com3 = log2(T3 / p2);
                    }
                }
            }
        }

        #pragma omp critical
        {
            if (found_local && (!found_any || best_local < best_cost)) {
                found_any = true;
                best_cost = best_local;
                best_global = bestL;
            }
        }
    }

    free(mu_of_eta);
    free(f_of_eta);
    free(p_of_eta);
    free(D_of_eta);
    free(m_of_eta);
    free(gamma_pre);
    free(removed_pre);

    if (!found_any) return 0;
    *out_best = best_global;
    return 1;
}

static void usage(const char *prog) {
    printf("Usage: %s n k h [binFlag=0 or 1] [xiCoef]\n", prog);
    printf("Example: %s 1024 652 106 0 1\n", prog);
#ifdef _OPENMP
    printf("OpenMP: enabled. Control threads with OMP_NUM_THREADS.\n");
#else
    printf("OpenMP: NOT enabled (recompile with -fopenmp).\n");
#endif
}

int main(int argc, char **argv) {
    if (argc < 4 || argc > 6) {
        usage(argv[0]);
        return 1;
    }

    int n = atoi(argv[1]);
    int k = atoi(argv[2]);
    int h = atoi(argv[3]);

    int binFlag = (argc >= 5) ? atoi(argv[4]) : 0;
    double xiCoef = (argc >= 6) ? atof(argv[5]) : (binFlag == 0)? 1.0 : 0.5;

    if (n <= 0 || k <= 0 || h <= 0) {
        printf("Error: require n,k,h > 0\n");
        return 1;
    }

    BestResultSearchParams best = {0};
    int ok = search_params_best_omp(n, k, h, 2.8, binFlag, xiCoef, &best);

    if (!ok) {
        printf("No candidate found.\n");
        return 0;
    }

    printf("Best log2(cost)=%.10f\n", best.best_log2_cost);
    printf("\n");

    printf("Truncated inputs)\n");
    printf("n=%d k=%d h=%d (numTrun=%d)\n", best.n_adj, best.k_adj, best.h, best.numTrun);
    printf("beta=%d omega=%d flin=%d\n", best.beta, best.omega, best.flin);
    printf("\n");

    printf("Best algorithmic parameters)\n");
    printf("(eta1,eta,eta_pre)\t\t=(%d,%d,%d)\n", best.eta1, best.eta1+best.eta2, best.eta_pre);
    printf("(f1,mu1,f,mu,f_pre,nu)\t=(%d,%d,%d,%d,%d,%d)\n", best.f1, best.mu1, best.f2, best.mu2, best.f_pre, best.nu);
    printf("\n");
    
    printf("Related values)\n");
    printf("n1=%d, m1=%lld, N1=%lld, D1=%lld, d1=%lld\n", best.n1, (long long)best.m1, (long long)best.N1, (long long)best.D1, (long long)best.d1);
    printf("n2=%d, m2=%lld, N2=%lld, D2=%lld, d2=%lld\n", best.n2, (long long)best.m2, (long long)best.N2, (long long)best.D2, (long long)best.d2);
    printf("gamma=%lld\n", (long long)best.gamma);
    printf("n2-eta_pre = %d\n", best.n2 - best.eta_pre);
    printf("\n");

    printf("Restriction values) restriction: N2pp x q^{-xi(N2pp)} << p_{eta}\n");
    printf("p_{eta}\t\t\t=2^{%.2f}\n", log2(best.p2));
    if (binFlag) {
        //printf("N2pp=%lld\n", (long long)best.N2pp);
        //printf("log N2pp=%.2f\n", log2((double)best.N2pp));
        //printf("log q^-xi(N2pp)=%.2f\n", (double)(-xi(best.N2pp, xiCoef, 2)));
        printf("N2pp x q^{-xi(N2pp)}\t=2^{%.2f}\n", log2((double)best.N2pp)-(double)xi(best.N2pp,xiCoef,2));
    }
    else {
        //printf("N2pp=%lld\n", (long long)best.N2pp);
        //printf("log N2pp=%.2f\n", log2((double)best.N2pp));
        //printf("log q^-xi(N2pp)=%.2f\n", (double)(-128*xi(best.N2pp, xiCoef, 3)));
        printf("N2pp x q^{-xi(N2pp)}\t=2^{%.2f}\n", log2((double)best.N2pp)-(double)(128*xi(best.N2pp,xiCoef,3)));
    }

    return 0;
}