#include "fft.h"

using namespace std;

template <int N, int S, int SO = 1, typename T=float>
struct prunefft {
  typedef complex<T> cfloat_t;
  prunefft<N/2, S*2, SO> s2;
  prunefft<N/4, S*4, SO> s4;
  void apply(const cfloat_t *in, cfloat_t *out, const cfloat_t *l1, const cfloat_t *l3, unsigned log2stride) {
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

template <int S, int SO, typename T>
struct prunefft<4, S, SO, T> {
  typedef complex<T> cfloat_t;
  void apply(const cfloat_t *in, cfloat_t *out, const cfloat_t *l1, const cfloat_t *l3, unsigned log2stride) {
    (void)l1;
    (void)l3;
    (void)log2stride;
    cfloat_t Uk = in[0];
    out[0] = Uk;
    out[SO*2] = Uk;
    out[SO] = Uk;
    out[SO*3] = Uk;
  }
};

template <int S, int SO, typename T>
struct prunefft<8, S, SO, T> {
  typedef complex<T> cfloat_t;
  void apply(const cfloat_t *in, cfloat_t *out, const cfloat_t *l1, const cfloat_t *l3, unsigned log2stride) {
    (void)l1;
    (void)l3;
    (void)log2stride;
    {
      cfloat_t Uk = in[0];
      out[SO*0] = Uk;
      out[SO*1] = Uk;
      out[SO*2] = Uk;
      out[SO*3] = Uk;
      out[SO*4] = Uk;
      out[SO*5] = Uk;
      out[SO*6] = Uk;
      out[SO*7] = Uk;
    }
  }
};
template <int S, int SO, typename T>
struct prunefft<16, S, SO, T> {
  typedef complex<T> cfloat_t;
  void apply(const cfloat_t *in, cfloat_t *out, const cfloat_t *l1, const cfloat_t *l3, unsigned log2stride) {
    (void)l1;
    (void)l3;
    (void)log2stride;
    {
      cfloat_t Uk = in[0];
      out[SO*0] = Uk;
      out[SO*1] = Uk;
      out[SO*2] = Uk;
      out[SO*3] = Uk;
      out[SO*4] = Uk;
      out[SO*5] = Uk;
      out[SO*6] = Uk;
      out[SO*7] = Uk;
      out[SO*8] = Uk;
      out[SO*9] = Uk;
      out[SO*10] = Uk;
      out[SO*11] = Uk;
      out[SO*12] = Uk;
      out[SO*13] = Uk;
      out[SO*14] = Uk;
      out[SO*15] = Uk;
    }
  }
};
template <int S, int SO, typename T>
struct prunefft<32, S, SO, T> {
  typedef complex<T> cfloat_t;
  void apply(const cfloat_t *in, cfloat_t *out, const cfloat_t *l1, const cfloat_t *l3, unsigned log2stride) {
    (void)l1;
    (void)l3;
    (void)log2stride;
    {
      cfloat_t Uk = in[0];
      out[SO*0] = Uk;
      out[SO*1] = Uk;
      out[SO*2] = Uk;
      out[SO*3] = Uk;
      out[SO*4] = Uk;
      out[SO*5] = Uk;
      out[SO*6] = Uk;
      out[SO*7] = Uk;
      out[SO*8] = Uk;
      out[SO*9] = Uk;
      out[SO*10] = Uk;
      out[SO*11] = Uk;
      out[SO*12] = Uk;
      out[SO*13] = Uk;
      out[SO*14] = Uk;
      out[SO*15] = Uk;
      out[SO*16] = Uk;
      out[SO*17] = Uk;
      out[SO*18] = Uk;
      out[SO*19] = Uk;
      out[SO*20] = Uk;
      out[SO*21] = Uk;
      out[SO*22] = Uk;
      out[SO*23] = Uk;
      out[SO*24] = Uk;
      out[SO*25] = Uk;
      out[SO*26] = Uk;
      out[SO*27] = Uk;
      out[SO*28] = Uk;
      out[SO*29] = Uk;
      out[SO*30] = Uk;
      out[SO*31] = Uk;
    }
  }
};

int main()
{
  complex<float> m[1024];
  complex<float> out[1024];
  complex<float> out2[1024];
  for (int i = 0; i != 16; ++i) {
    //m[i] = {float(i+1000), float(2*i+100)};
    m[1023-i] = float(i+2000);
  }
  const lut<256, 1> l1;
  const lut<256, 3> l3;
  splitfft<1024, 1, 1> sfft;
  prunefft<1024, 1, 1> pfft;
  sfft.apply(m, out, l1.d, l3.d, 0);
  pfft.apply(m, out2, l1.d, l3.d, 0);
  //printf("%f,%f\n", static_cast<double>(cshift), static_cast<double>(rshift));
  for (int i = 0; i != 1024; ++i) {
    if (abs(out[i]-out2[i]) > 1e-7)
      printf("prune failed at %d!\n", i);
  }
  return 0;
}
