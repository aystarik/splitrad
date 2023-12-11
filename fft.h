#pragma once

#include <cmath>
#include <complex>

template <int N, int k, typename T = float>
struct lut {
    typedef std::complex<T> cfloat_t;
    cfloat_t d[N];
    lut() {
        constexpr cfloat_t w0 = {0.0, -T(M_PI)/T(2 * N)};
        d[0] = {1.0, 0.0};
        for (int i = 1; i < N; ++i)
            d[i] = exp(w0 * float(i * k));
    }
};

// generic N = 2^a case
template <int N, int S, int SO = 1, typename T = float>
struct splitfft {
    typedef std::complex<T> cfloat_t;
    splitfft<N / 2, S * 2, SO> s2;
    splitfft<N / 4, S * 4, SO> s4;
    void apply(const cfloat_t * __restrict__ in, cfloat_t *out,
               const cfloat_t * __restrict__ l1, const cfloat_t * __restrict__ l3,
               unsigned log2stride) {
        constexpr cfloat_t I(0.0f, 1.0f);
        s2.apply(in, out, l1, l3, log2stride + 1);
        s4.apply(in + S, out + N/2*SO, l1, l3, log2stride + 2);
        s4.apply(in + 3*S, out + 3*N/4*SO, l1, l3, log2stride + 2);
        {
            cfloat_t Uk = out[0];
            cfloat_t Zk = out[0+SO*N/2];
            cfloat_t Uk2 = out[0+SO*N/4];
            cfloat_t Zdk = out[0+SO*3*N/4];
            out[0] = Uk + (Zk + Zdk);
            out[0+SO*N/2] = Uk - (Zk + Zdk);
            out[0+SO*N/4] = Uk2 - I*(Zk - Zdk);
            out[0+SO*3*N/4] = Uk2 + I*(Zk - Zdk);
        }

        for(unsigned k = 1; k < N/4; ++k) {
            cfloat_t Uk = out[SO*(k)];
            cfloat_t Zk = out[SO*(k+N/2)];
            cfloat_t Uk2 = out[SO*(k+N/4)];
            cfloat_t Zdk = out[SO*(k+3*N/4)];
            cfloat_t w1 = l1[k << log2stride];
            cfloat_t w3 = l3[k << log2stride];
            out[SO*(k)] = Uk + (w1*Zk + w3*Zdk);
            out[SO*(k+N/2)] = Uk - (w1*Zk + w3*Zdk);
            out[SO*(k+N/4)] = Uk2 - I*(w1*Zk - w3*Zdk);
            out[SO*(k+3*N/4)] = Uk2 + I*(w1*Zk - w3*Zdk);
        }
    }
};

// template specialization for N = 4
template <int S, int SO, typename T>
struct splitfft<4, S, SO, T> {
    typedef std::complex<T> cfloat_t;
    void apply(const cfloat_t * __restrict__ in, cfloat_t *out,
               const cfloat_t * __restrict__ l1, const cfloat_t * __restrict__ l3,
               unsigned log2stride) {
        (void)l1;
        (void)l3;
        (void)log2stride;
        constexpr cfloat_t I(0.0f, 1.0f);

        cfloat_t a =   in[0];
        cfloat_t Zk =  in[S];
        cfloat_t b =   in[2 * S];
        cfloat_t Zdk = in[3 * S];
        cfloat_t Uk = a + b;
        cfloat_t Uk2 = a - b;

        out[0] = Uk + (Zk + Zdk);
        out[SO] = Uk2 - I * (Zk - Zdk);
        out[SO * 2] = Uk - (Zk + Zdk);
        out[SO * 3] = Uk2 + I * (Zk - Zdk);
    }
};

// template specialization for N = 8
template <int S, int SO, typename T>
struct splitfft<8, S, SO, T> {
    typedef std::complex<T> cfloat_t;
    void apply(const cfloat_t * __restrict__ in, cfloat_t *out,
               const cfloat_t * __restrict__ l1, const cfloat_t * __restrict__ l3, unsigned log2stride) {
        (void)l1;
        (void)l3;
        (void)log2stride;
        constexpr cfloat_t I(0.0f, 1.0f);
        {
            cfloat_t a =   in[0];
            cfloat_t Zk =  in[2 * S];
            cfloat_t b =   in[4 * S];
            cfloat_t Zdk = in[6 * S];
            cfloat_t Uk = a + b;
            cfloat_t Uk2 = a - b;
            a = Zk + Zdk;
            b = I * (Zk - Zdk);
            out[0] = Uk + a;
            out[SO] = Uk2 - b;
            out[SO * 2] = Uk - a;
            out[SO * 3] = Uk2 + b;
            a =          in[S];
            cfloat_t c = in[3 * S];
            b =          in[5 * S];
            cfloat_t d = in[7 * S];
            out[4 * SO] = a + b;
            out[5 * SO] = a - b;
            out[6 * SO] = c + d;
            out[7 *SO] = c - d;
        }
        {
            cfloat_t Uk =  out[0];
            cfloat_t Uk2 = out[SO * 2];
            cfloat_t Zk =  out[SO * 4];
            cfloat_t Zdk = out[SO * 6];
            out[0] = Uk + (Zk + Zdk);
            out[SO * 2] = Uk2 - I*(Zk - Zdk);
            out[SO * 4] = Uk - (Zk + Zdk);
            out[SO * 6] = Uk2 + I*(Zk - Zdk);
        }
        {
            cfloat_t Uk =  out[SO * 1];
            cfloat_t Uk2 = out[SO * 3];
            cfloat_t Zk =  out[SO * 5];
            cfloat_t Zdk = out[SO * 7];
            constexpr float w = M_SQRT1_2;
            constexpr cfloat_t w1 = {w, -w};
            constexpr cfloat_t w3 = {-w, -w};
            cfloat_t a = (w1 * Zk + w3 * Zdk);
            cfloat_t b = I * (w1 * Zk - w3 * Zdk);
            out[SO * 1] = Uk + a;
            out[SO * 3] = Uk2 - b;
            out[SO * 5] = Uk - a;
            out[SO * 7] = Uk2 + b;
        }
    }
};

// 1D recursive template split-radix out-of-place FFT
template <int N, int S = 1, typename T = float>
struct fft1d {
    const lut<N/4, 1> l1;
    const lut<N/4, 3> l3;
    splitfft<N, S, S> sfft;
    void apply(const std::complex<T> *__restrict__ in, std::complex<T> *out) {
        sfft.apply(in, out, l1.d, l3.d, 0);
    }
};


// 2D recursive template split-radix out-of-place FFT
template <int N, typename T = float>
struct fft2d {
    typedef std::complex<T> cfloat_t;
    fft1d<N, N> fftc;
    fft1d<N, 1> fftr;
    cfloat_t temp[N*N];
    void apply(cfloat_t *in) {
        int c = N;
        cfloat_t *t = in;
        cfloat_t *o = temp;
        while (c--)
            fftc.apply(t++, o++);
        c = N;
        t = temp;
        o = in;
        while (c--) {
            fftr.apply(t, o);
            t += N;
            o += N;
        }
    }
    void apply(const cfloat_t *__restrict__ in, cfloat_t *out) {
        const cfloat_t *t = in;
        cfloat_t *o = temp;
        int c = N;
        while (c--)
            fftc.apply(t++, o++);
        t = temp;
        o = out;
        c = N;
        while (c--) {
            fftr.apply(t, o);
            t += N;
            o += N;
        }
    }
};
