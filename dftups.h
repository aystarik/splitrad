#pragma once

#include "fft.h"

template <unsigned N, typename T=float>
void conjmult(const std::complex<T> *a, const std::complex<T> *b, std::complex<T> *c) {
    const T *ai = reinterpret_cast<const T*>(a);
    const T *bi = reinterpret_cast<const T*>(b);
    T *ci = reinterpret_cast<T*>(c);
    int i = N * N;
    while (i--) {
        T air = ai[0];
        T aii = ai[1];
        T bir = bi[0];
        T bii = bi[1];
        ci[0] = air * bir + aii * bii;
        ci[1] = -aii * bir + air * bii;
        ci += 2, ai += 2, bi += 2;
    }
}

template <unsigned N, typename T=float>
void absmax(const std::complex<T> *m, int &c, int &r) {
    T maxv = -1.0f;
    unsigned maxi = 0;
    const std::complex<T> *mi = m;
    for (unsigned i = 0; i < N*N; ++i) {
        T x = mi[0].real();
        T y = mi[0].imag();
        mi++;
        T cabs = x*x + y*y;
        if (cabs > maxv) {
            maxv = cabs;
            maxi = i;
        }
    }
    c = maxi/N;
    r = maxi & (N-1);
}

template <int N, typename T=float>
struct twiddle {
    T tw[N*2];
    twiddle() {
        constexpr float k = M_PIf32/static_cast<T>(N);
        tw[N/2 + 1] = tw[0] = 1.0f;
        tw[N/2] = tw[3*N/2] = tw[N + 1] = tw[1] = 0.0f;
        tw[N] = tw[3*N/2 + 1] = -1.0f;
        for (unsigned i = 2; i < N/2; i += 2) {
            T c = cos(k*i);
            tw[N/2 - i + 1] = tw[N/2 + i + 1] = c;
            tw[3*N/2 - i + 1] = tw[3*N/2 + i + 1] = -c;
            tw[i] = tw[2*N - i] = c;
            tw[N - i] = tw[N + i] = -c;
        }
    }
    std::complex<T> sincos(int32_t i) const {
        int32_t mask = i >> 31;
        const T *p = tw + (((N-1)&((mask + i)^(mask))) << 1);
        T c = p[0];
        T x = p[1];
        T s = (mask)?x:-x;
        return {c,s};
    }
};

template <int N, typename T=float>
void matmult(const std::complex<T> *a, const std::complex<T> *b, std::complex<T> *c) {
    std::complex<T> *cij = c;
    const std::complex<T> *ai = a;
    int i = N;
    while (i--) {
        const std::complex<T> *bj = b;
        int j = N;
        while (j--) {
            std::complex<T> t{0.0, 0.0};
            const std::complex<T> *aik = ai, *bkj = bj;
            int k = N/2;
            while (k--) {
                t += aik[0] * bkj[0*N];
                t += aik[1] * bkj[1*N];
                aik += 2;
                bkj += 2*N;
            }
            ++bj;
            *cij++ = t;
        }
        ai += N;
    }
}

template <int N, typename T=float>
void matmultmax(const std::complex<T> *a, const std::complex<T> *b, int &coff, int &roff) {
    const std::complex<T> *ai = a;
    T maxv = -1.0f;
    for (int i = 0; i < N; ++i) {
        const std::complex<T> *bj = b;
        for (int j = 0; j < N; ++j) {
            std::complex<T> t{0.0, 0.0};
            const std::complex<T> *aik = ai, *bkj = bj;
            int k = N/2;
            while (k--) {
                t += aik[0] * bkj[0*N];
                t += aik[1] * bkj[1*N];
                aik += 2;
                bkj += 2*N;
            }
            T x = t.real();
            T y = t.imag();
            T cabs = x*x + y*y;
            if (cabs > maxv) {
                maxv = cabs;
                coff = i;
                roff = j;
            }
            ++bj;
        }
        ai += N;
    }
}

template <int N, int K, typename T = float>
struct dftups {
    typedef twiddle<N*K, T> TW;
    const TW tw;
    void apply(const std::complex<T> *in, int &coff, int &roff) {
        std::complex<T> kern[N*N], temp[N*N];
        std::complex<T> *k1 = kern;
        int cnt1 = N;
        int step = -coff;
        while(cnt1--) {
            int a1 = 0, a2 = -N/2 * step;
            int cnt2 = N/2;
            while (cnt2--) {
                k1[0] = tw.sincos(a1);
                a1 += step;
                k1[N/2] = tw.sincos(a2);
                ++k1;
                a2 += step;
            }
            k1 += N/2;
            ++step;
        }
        matmult<N,T>(kern, in, temp);
        k1 = kern;
        int step1 = 0;
        int step2 = -N/2;
        int i = N/2;
        while (i--) {
            int a1 = -step1 * roff;
            int a2 = -step2 * roff;
            int j = N;
            while (j--) {
                k1[0] = tw.sincos(a1);
                a1 += step1;
                k1[N*N/2] = tw.sincos(a2);
                ++k1;
                a2 += step2;
            }
            ++step1;
            ++step2;
        }
        matmultmax<N,T>(temp, kern, coff, roff);
    }
};
