#include "fft.h"
#include "dftups.h"

#define NEED_TEST_DATA
#include "test_data.cc"

template <unsigned N, typename T=float>

void im2float(const uint16_t *d, std::complex<T> *tt) {
    constexpr float k = 1.0f/65536.0f;
    int i = N*N;
    while (i--) {
        *tt++ = {k * *d++, 0.0};
    }
}

void test(float &c, float &r)
{
    constexpr int N = 32, NN = N*N, K = 16;
    fft2d<N> fft;
    typedef std::complex<float> cfloat_t;
    //float rshift = 0, cshift = 0;
    cfloat_t im1[NN], im2[NN];
    im2float<N>(b, im1);
    fft.apply(im1);
    dftups<N, K> dups;
    cfloat_t cc[NN], CC[NN];
    for (unsigned i = 0; i < 1; ++i) {
        //palToggleLine(LINE_LED1);
        im2float<N>(a, im2);
        fft.apply(im2);
        conjmult<N>(im1, im2, cc);
        fft.apply(cc, CC);
        int coff = 0, roff = 0;
        absmax<N>(CC, coff, roff);
        if (coff > N/2)
            coff -= N;
        if (roff > N/2)
            roff -= N;
        int coffu = N/2 - coff*K;
        int roffu = N/2 - roff*K;
        dups.apply(cc, coffu, roffu);
        roffu -= N/2;
        coffu -= N/2;
        r = roff + roffu/float(K);
        c = coff + coffu/float(K);
    }
}

int main()
{
  float cshift, rshift;
  test(cshift, rshift);
  printf("%f,%f\n", static_cast<double>(cshift), static_cast<double>(rshift));
  return 0;
}
