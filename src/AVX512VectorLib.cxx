#include "AVX512VectorLib.hpp"



namespace InstabilityDetectionLib {
//typedef AVL_IntegerVectorType AVL_IntegerVectorType; // vector of 16 int   (avx)

/* yes I know, the top of this file is quite ugly */
#ifdef _MSC_VER /* visual c++ */
# define ALIGN32_BEG __declspec(align(32))
# define ALIGN32_END
#else /* gcc or icc */
# define ALIGN32_BEG
# define ALIGN32_END __attribute__((aligned(32)))
#endif


_PS512_CONST(1  , 1.0f);
_PS512_CONST(2  , 2.0f);
_PS512_CONST(0p5, 0.5f);
/* the smallest non denormalized float number */
_PS512_CONST_TYPE(min_norm_pos, int, 0x00800000);
//_PS512_CONST_TYPE(mant_mask, int, 0x7f800000);
_PS512_CONST_TYPE(inv_mant_mask, int, ~0x7f800000);

_PS512_CONST_TYPE(sign_mask, int, (int)0x80000000);
_PS512_CONST_TYPE(inv_sign_mask, int, ~0x80000000);

_PI32_CONST512(0, 0);
_PI32_CONST512(1, 1);
_PI32_CONST512(0xffffffff, (int)0xFFFFFFFF);
_PI32_CONST512(inv1, ~1);
_PI32_CONST512(2, 2);
_PI32_CONST512(4, 4);
_PI32_CONST512(0x7f, 0x7f);

_PS512_CONST(cephes_SQRTHF, 0.707106781186547524f);
_PS512_CONST(cephes_log_p0, 7.0376836292E-2f);
_PS512_CONST(cephes_log_p1, - 1.1514610310E-1f);
_PS512_CONST(cephes_log_p2, 1.1676998740E-1f);
_PS512_CONST(cephes_log_p3, - 1.2420140846E-1f);
_PS512_CONST(cephes_log_p4, + 1.4249322787E-1f);
_PS512_CONST(cephes_log_p5, - 1.6668057665E-1f);
_PS512_CONST(cephes_log_p6, + 2.0000714765E-1f);
_PS512_CONST(cephes_log_p7, - 2.4999993993E-1f);
_PS512_CONST(cephes_log_p8, + 3.3333331174E-1f);
_PS512_CONST(cephes_log_q1, -2.12194440e-4f);
_PS512_CONST(cephes_log_q2, 0.693359375f);

//static inline AVL_IntegerVectorType _wrapAVL_slli_epi32(AVL_IntegerVectorType x, int  y) { return AVL_slli_epi32(x, y); }
//static inline AVL_IntegerVectorType _wrapAVL_srli_epi32(AVL_IntegerVectorType x, int  y) { return AVL_srli_epi32(x, y); }
//static inline AVL_IntegerVectorType _wrapAVL_sub_epi32 (AVL_IntegerVectorType x, AVL_IntegerVectorType y) { return AVL_sub_epi32 (x, y); }
//static inline AVL_IntegerVectorType _wrapAVL_add_epi32 (AVL_IntegerVectorType x, AVL_IntegerVectorType y) { return AVL_add_epi32 (x, y); }


/* natural logarithm computed for 16 simultaneous float
   return NaN for x <= 0
*/

AVL_FloatVectorType log512_ps(AVL_FloatVectorType x) {
  AVL_IntegerVectorType imm0;
  AVL_FloatVectorType one = *(AVL_FloatVectorType*)_ps512_1;

  x = AVL_max_ps(x, *(AVL_FloatVectorType*)_ps512_min_norm_pos);  /* cut off denormalized stuff */

  // can be done with AVX2
  imm0 = AVL_srli_epi32(AVL_castps_si(x), 23);

  /* keep only the fractional part */
  x = AVL_castsi_ps(AVL_and_si(AVL_castps_si(x), AVL_castps_si(*(AVL_FloatVectorType*)_ps512_inv_mant_mask)));
  x = AVL_castsi_ps(AVL_or_si(AVL_castps_si(x), AVL_castps_si(*(AVL_FloatVectorType*)_ps512_0p5)));

  // this is again another AVX2 instruction
  imm0 = AVL_sub_epi32(imm0, *(AVL_IntegerVectorType*)_pi32_512_0x7f);
  AVL_FloatVectorType e = AVL_cvtepi32_ps(imm0);

  e = AVL_add_ps(e, one);

  /* part2:
     if( x < SQRTHF ) {
       e -= 1;
       x = x + x - 1.0;
     } else { x = x - 1.0; }
  */
  AVL_MaskType mask2 = AVL_cmp_ps_mask(x, *(AVL_FloatVectorType*)_ps512_cephes_SQRTHF, _CMP_LT_OS);
  e = AVL_mask_sub_ps(e, mask2, e, one);
  x = AVL_mask_fmsub_ps(x, mask2, *(AVL_FloatVectorType*)_ps512_2, *(AVL_FloatVectorType*)_ps512_1);
  mask2 = AVL_knot(mask2);
  x = AVL_mask_sub_ps(x, mask2, x, one);


  AVL_FloatVectorType z = AVL_mul_ps(x, x);

  AVL_FloatVectorType y = AVL_fmadd_ps(*(AVL_FloatVectorType*)_ps512_cephes_log_p0, x, *(AVL_FloatVectorType*)_ps512_cephes_log_p1); //0
  y = AVL_fmadd_ps(y, x, *(AVL_FloatVectorType*)_ps512_cephes_log_p2); //1
  y = AVL_fmadd_ps(y, x, *(AVL_FloatVectorType*)_ps512_cephes_log_p3); //2
  y = AVL_fmadd_ps(y, x, *(AVL_FloatVectorType*)_ps512_cephes_log_p4); //3
  y = AVL_fmadd_ps(y, x, *(AVL_FloatVectorType*)_ps512_cephes_log_p5); //4
  y = AVL_fmadd_ps(y, x, *(AVL_FloatVectorType*)_ps512_cephes_log_p6); //5
  y = AVL_fmadd_ps(y, x, *(AVL_FloatVectorType*)_ps512_cephes_log_p7); //6
  y = AVL_fmadd_ps(y, x, *(AVL_FloatVectorType*)_ps512_cephes_log_p8); //7

  y = AVL_mul_ps(y, x);
  y = AVL_mul_ps(y, z);

  y = AVL_fmadd_ps(e, *(AVL_FloatVectorType*)_ps512_cephes_log_q1, y);

  y = AVL_fnmadd_ps(z, *(AVL_FloatVectorType*)_ps512_0p5, y);

  x = AVL_add_ps(x, y);

  x = AVL_fmadd_ps(e, *(AVL_FloatVectorType*)_ps512_cephes_log_q2, x);

  return x;
}




_PS512_CONST(exp_hi,  88.3762626647950f);

_PS512_CONST(exp_lo,  -88.3762626647949f);

_PS512_CONST(cephes_LOG2EF, 1.44269504088896341f);
_PS512_CONST(cephes_exp_C1, 0.693359375f);
_PS512_CONST(cephes_exp_C2, -2.12194440e-4f);

_PS512_CONST(cephes_exp_p0, 1.9875691500E-4f);
_PS512_CONST(cephes_exp_p1, 1.3981999507E-3f);
_PS512_CONST(cephes_exp_p2, 8.3334519073E-3f);
_PS512_CONST(cephes_exp_p3, 4.1665795894E-2f);
_PS512_CONST(cephes_exp_p4, 1.6666665459E-1f);
_PS512_CONST(cephes_exp_p5, 5.0000001201E-1f);

AVL_FloatVectorType exp512_ps(AVL_FloatVectorType x) {
  AVL_FloatVectorType tmp = AVL_setzero_ps(), fx;
  AVL_IntegerVectorType imm0;
  AVL_FloatVectorType one = *(AVL_FloatVectorType*)_ps512_1;

  x = AVL_min_ps(x, *(AVL_FloatVectorType*)_ps512_exp_hi);
  x = AVL_max_ps(x, *(AVL_FloatVectorType*)_ps512_exp_lo);

  /* express exp(x) as exp(g + n*log(2)) */
  fx = AVL_mul_ps(x, *(AVL_FloatVectorType*)_ps512_cephes_LOG2EF);
  fx = AVL_add_ps(fx, *(AVL_FloatVectorType*)_ps512_0p5);

  /* how to perform a floorf with SSE: just below */
  //imm0 = AVL_cvttps_epi32(fx);
  //tmp  = AVL_cvtepi32_ps(imm0);

  tmp = AVL_floor_ps(fx);

  /* if greater, substract 1 */
  AVL_MaskType mask2 = AVL_cmp_ps_mask(tmp, fx, _CMP_GT_OS);
  AVL_FloatVectorType mask = AVL_mask_blend_ps(mask2, *(AVL_FloatVectorType*)_pi32_512_0, *(AVL_FloatVectorType*)_pi32_512_0xffffffff);
  mask = AVL_castsi_ps(AVL_and_si(AVL_castps_si(mask), AVL_castps_si(one)));
  fx = AVL_sub_ps(tmp, mask);

  tmp = AVL_mul_ps(fx, *(AVL_FloatVectorType*)_ps512_cephes_exp_C1);
  AVL_FloatVectorType z = AVL_mul_ps(fx, *(AVL_FloatVectorType*)_ps512_cephes_exp_C2);
  x = AVL_sub_ps(x, tmp);
  x = AVL_sub_ps(x, z);

  z = AVL_mul_ps(x, x);

  AVL_FloatVectorType y = *(AVL_FloatVectorType*)_ps512_cephes_exp_p0;
  y = AVL_mul_ps(y, x);
  y = AVL_add_ps(y, *(AVL_FloatVectorType*)_ps512_cephes_exp_p1);
  y = AVL_mul_ps(y, x);
  y = AVL_add_ps(y, *(AVL_FloatVectorType*)_ps512_cephes_exp_p2);
  y = AVL_mul_ps(y, x);
  y = AVL_add_ps(y, *(AVL_FloatVectorType*)_ps512_cephes_exp_p3);
  y = AVL_mul_ps(y, x);
  y = AVL_add_ps(y, *(AVL_FloatVectorType*)_ps512_cephes_exp_p4);
  y = AVL_mul_ps(y, x);
  y = AVL_add_ps(y, *(AVL_FloatVectorType*)_ps512_cephes_exp_p5);
  y = AVL_mul_ps(y, z);
  y = AVL_add_ps(y, x);
  y = AVL_add_ps(y, one);

  /* build 2^n */
  imm0 = AVL_cvttps_epi32(fx);
  // another two AVX2 instructions
  imm0 = AVL_add_epi32(imm0, *(AVL_IntegerVectorType*)_pi32_512_0x7f);
  imm0 = AVL_slli_epi32(imm0, 23);
  AVL_FloatVectorType pow2n = AVL_castsi_ps(imm0);
  y = AVL_mul_ps(y, pow2n);
  return y;
}

_PS512_CONST(minus_cephes_DP1, -0.78515625f);
_PS512_CONST(minus_cephes_DP2, -2.4187564849853515625e-4f);
_PS512_CONST(minus_cephes_DP3, -3.77489497744594108e-8f);
_PS512_CONST(sincof_p0, -1.9515295891E-4f);
_PS512_CONST(sincof_p1,  8.3321608736E-3f);
_PS512_CONST(sincof_p2, -1.6666654611E-1f);
_PS512_CONST(coscof_p0,  2.443315711809948E-005f);
_PS512_CONST(coscof_p1, -1.388731625493765E-003f);
_PS512_CONST(coscof_p2,  4.166664568298827E-002f);
_PS512_CONST(cephes_FOPI, 1.27323954473516f); // 4 / M_PI


/* evaluation of 16 sines at onces using AVX intrisics

   The code is the exact rewriting of the cephes sinf function.
   Precision is excellent as long as x < 8192 (I did not bother to
   take into account the special handling they have for greater values
   -- it does not return garbage for arguments over 8192, though, but
   the extra precision is missing).

   Note that it is such that sinf((float)M_PI) = 8.74e-8, which is the
   surprising but correct result.

*/
AVL_FloatVectorType sin512_ps(AVL_FloatVectorType x) { // any x
  AVL_FloatVectorType xmm1, xmm2 = AVL_setzero_ps(), xmm3, sign_bit, y;
  AVL_IntegerVectorType imm0, imm2;

  sign_bit = x;
  /* take the absolute value */

  x = AVL_castsi_ps(AVL_and_si(AVL_castps_si(x), AVL_castps_si(*(AVL_FloatVectorType*)_ps512_inv_sign_mask)));
  /* extract the sign bit (upper one) */

  sign_bit = AVL_castsi_ps(AVL_and_si(AVL_castps_si(sign_bit), AVL_castps_si(*(AVL_FloatVectorType*)_ps512_sign_mask)));

  /* scale by 4/Pi */
  y = AVL_mul_ps(x, *(AVL_FloatVectorType*)_ps512_cephes_FOPI);

  /*
    Here we start a series of integer operations, which are in the
    realm of AVX2.
    If we don't have AVX, let's perform them using SSE2 directives
  */

  /* store the integer part of y in mm0 */
  imm2 = AVL_cvttps_epi32(y);
  /* j=(j+1) & (~1) (see the cephes sources) */
  // another two AVX2 instruction
  imm2 = AVL_add_epi32(imm2, *(AVL_IntegerVectorType*)_pi32_512_1);
  imm2 = AVL_and_si(imm2, *(AVL_IntegerVectorType*)_pi32_512_inv1);
  y = AVL_cvtepi32_ps(imm2);

  /* get the swap sign flag */
  imm0 = AVL_and_si(imm2, *(AVL_IntegerVectorType*)_pi32_512_4);
  imm0 = AVL_slli_epi32(imm0, 29);
  /* get the polynom selection mask
     there is one polynom for 0 <= x <= Pi/4
     and another one for Pi/4<x<=Pi/2

     Both branches will be computed.
  */
  imm2 = AVL_and_si(imm2, *(AVL_IntegerVectorType*)_pi32_512_2);


  AVL_MaskType imm22 = AVL_cmpeq_epi32_mask(imm2, *(AVL_IntegerVectorType*)_pi32_512_0);
  imm2 = AVL_mask_blend_epi32(imm22, *(AVL_IntegerVectorType*)_pi32_512_0, *(AVL_IntegerVectorType*)_pi32_512_0xffffffff);

  AVL_FloatVectorType swap_sign_bit = AVL_castsi_ps(imm0);
  AVL_FloatVectorType poly_mask = AVL_castsi_ps(imm2);

  sign_bit = AVL_castsi_ps(AVL_xor_si(AVL_castps_si(sign_bit), AVL_castps_si(swap_sign_bit)));

  /* The magic pass: "Extended precision modular arithmetic"
     x = ((x - y * DP1) - y * DP2) - y * DP3; */
  xmm1 = *(AVL_FloatVectorType*)_ps512_minus_cephes_DP1;
  xmm2 = *(AVL_FloatVectorType*)_ps512_minus_cephes_DP2;
  xmm3 = *(AVL_FloatVectorType*)_ps512_minus_cephes_DP3;
  xmm1 = AVL_mul_ps(y, xmm1);
  xmm2 = AVL_mul_ps(y, xmm2);
  xmm3 = AVL_mul_ps(y, xmm3);
  x = AVL_add_ps(x, xmm1);
  x = AVL_add_ps(x, xmm2);
  x = AVL_add_ps(x, xmm3);

  /* Evaluate the first polynom  (0 <= x <= Pi/4) */
  y = *(AVL_FloatVectorType*)_ps512_coscof_p0;
  AVL_FloatVectorType z = AVL_mul_ps(x, x);

  y = AVL_mul_ps(y, z);
  y = AVL_add_ps(y, *(AVL_FloatVectorType*)_ps512_coscof_p1);
  y = AVL_mul_ps(y, z);
  y = AVL_add_ps(y, *(AVL_FloatVectorType*)_ps512_coscof_p2);
  y = AVL_mul_ps(y, z);
  y = AVL_mul_ps(y, z);
  AVL_FloatVectorType tmp = AVL_mul_ps(z, *(AVL_FloatVectorType*)_ps512_0p5);
  y = AVL_sub_ps(y, tmp);
  y = AVL_add_ps(y, *(AVL_FloatVectorType*)_ps512_1);

  /* Evaluate the second polynom  (Pi/4 <= x <= 0) */

  AVL_FloatVectorType y2 = *(AVL_FloatVectorType*)_ps512_sincof_p0;
  y2 = AVL_mul_ps(y2, z);
  y2 = AVL_add_ps(y2, *(AVL_FloatVectorType*)_ps512_sincof_p1);
  y2 = AVL_mul_ps(y2, z);
  y2 = AVL_add_ps(y2, *(AVL_FloatVectorType*)_ps512_sincof_p2);
  y2 = AVL_mul_ps(y2, z);
  y2 = AVL_mul_ps(y2, x);
  y2 = AVL_add_ps(y2, x);

  /* select the correct result from the two polynoms */
  xmm3 = poly_mask;
  y2 = AVL_castsi_ps(AVL_and_si(AVL_castps_si(xmm3), AVL_castps_si(y2)));
  y = AVL_castsi_ps(AVL_andnot_si(AVL_castps_si(xmm3), AVL_castps_si(y)));
  y = AVL_add_ps(y, y2);
  /* update the sign */
  y = AVL_castsi_ps(AVL_xor_si(AVL_castps_si(y), AVL_castps_si(sign_bit)));

  return y;
}


/* almost the same as sin_ps */
AVL_FloatVectorType cos512_ps(AVL_FloatVectorType x) { // any x
  AVL_FloatVectorType xmm1, xmm2 = AVL_setzero_ps(), xmm3, y;
  AVL_IntegerVectorType imm0, imm2;

  /* take the absolute value */
  x = AVL_castsi_ps(AVL_and_si(AVL_castps_si(x), AVL_castps_si(*(AVL_FloatVectorType*)_ps512_inv_sign_mask)));

  /* scale by 4/Pi */
  y = AVL_mul_ps(x, *(AVL_FloatVectorType*)_ps512_cephes_FOPI);

  /* store the integer part of y in mm0 */
  imm2 = AVL_cvttps_epi32(y);
  /* j=(j+1) & (~1) (see the cephes sources) */
  imm2 = AVL_add_epi32(imm2, *(AVL_IntegerVectorType*)_pi32_512_1);
  imm2 = AVL_and_si(imm2, *(AVL_IntegerVectorType*)_pi32_512_inv1);
  y = AVL_cvtepi32_ps(imm2);
  imm2 = AVL_sub_epi32(imm2, *(AVL_IntegerVectorType*)_pi32_512_2);

  /* get the swap sign flag */
  imm0 = AVL_andnot_si(imm2, *(AVL_IntegerVectorType*)_pi32_512_4);
  imm0 = AVL_slli_epi32(imm0, 29);
  /* get the polynom selection mask */
  imm2 = AVL_and_si(imm2, *(AVL_IntegerVectorType*)_pi32_512_2);
  AVL_MaskType imm22 = AVL_cmpeq_epi32_mask(imm2, *(AVL_IntegerVectorType*)_pi32_512_0);
  imm2 = AVL_mask_blend_epi32(imm22, *(AVL_IntegerVectorType*)_pi32_512_0, *(AVL_IntegerVectorType*)_pi32_512_0xffffffff);

  AVL_FloatVectorType sign_bit = AVL_castsi_ps(imm0);
  AVL_FloatVectorType poly_mask = AVL_castsi_ps(imm2);

  /* The magic pass: "Extended precision modular arithmetic"
     x = ((x - y * DP1) - y * DP2) - y * DP3; */
  xmm1 = *(AVL_FloatVectorType*)_ps512_minus_cephes_DP1;
  xmm2 = *(AVL_FloatVectorType*)_ps512_minus_cephes_DP2;
  xmm3 = *(AVL_FloatVectorType*)_ps512_minus_cephes_DP3;
  xmm1 = AVL_mul_ps(y, xmm1);
  xmm2 = AVL_mul_ps(y, xmm2);
  xmm3 = AVL_mul_ps(y, xmm3);
  x = AVL_add_ps(x, xmm1);
  x = AVL_add_ps(x, xmm2);
  x = AVL_add_ps(x, xmm3);

  /* Evaluate the first polynom  (0 <= x <= Pi/4) */
  y = *(AVL_FloatVectorType*)_ps512_coscof_p0;
  AVL_FloatVectorType z = AVL_mul_ps(x, x);

  y = AVL_mul_ps(y, z);
  y = AVL_add_ps(y, *(AVL_FloatVectorType*)_ps512_coscof_p1);
  y = AVL_mul_ps(y, z);
  y = AVL_add_ps(y, *(AVL_FloatVectorType*)_ps512_coscof_p2);
  y = AVL_mul_ps(y, z);
  y = AVL_mul_ps(y, z);
  AVL_FloatVectorType tmp = AVL_mul_ps(z, *(AVL_FloatVectorType*)_ps512_0p5);
  y = AVL_sub_ps(y, tmp);
  y = AVL_add_ps(y, *(AVL_FloatVectorType*)_ps512_1);

  /* Evaluate the second polynom  (Pi/4 <= x <= 0) */

  AVL_FloatVectorType y2 = *(AVL_FloatVectorType*)_ps512_sincof_p0;
  y2 = AVL_mul_ps(y2, z);
  y2 = AVL_add_ps(y2, *(AVL_FloatVectorType*)_ps512_sincof_p1);
  y2 = AVL_mul_ps(y2, z);
  y2 = AVL_add_ps(y2, *(AVL_FloatVectorType*)_ps512_sincof_p2);
  y2 = AVL_mul_ps(y2, z);
  y2 = AVL_mul_ps(y2, x);
  y2 = AVL_add_ps(y2, x);

  /* select the correct result from the two polynoms */
  xmm3 = poly_mask;
  y2 = AVL_castsi_ps(AVL_and_si(AVL_castps_si(xmm3), AVL_castps_si(y2)));
  y = AVL_castsi_ps(AVL_andnot_si(AVL_castps_si(xmm3), AVL_castps_si(y)));
  y = AVL_add_ps(y, y2);
  /* update the sign */
  y = AVL_castsi_ps(AVL_xor_si(AVL_castps_si(y), AVL_castps_si(sign_bit)));

  return y;
}

/* since sin512_ps and cos512_ps are almost identical, sincos512_ps could replace both of them..
   it is almost as fast, and gives you a free cosine with your sine */
void sincos512_ps(AVL_FloatVectorType x, AVL_FloatVectorType *s, AVL_FloatVectorType *c) {

  AVL_FloatVectorType xmm1, xmm2, xmm3 = AVL_setzero_ps(), sign_bit_sin, y;
  AVL_IntegerVectorType imm0, imm2, imm4;

  sign_bit_sin = x;
  /* take the absolute value */
  x = AVL_castsi_ps(AVL_and_si(AVL_castps_si(x), AVL_castps_si(*(AVL_FloatVectorType*)_ps512_inv_sign_mask)));
  /* extract the sign bit (upper one) */
  sign_bit_sin = AVL_castsi_ps(AVL_and_si(AVL_castps_si(sign_bit_sin), AVL_castps_si(*(AVL_FloatVectorType*)_ps512_sign_mask)));

  /* scale by 4/Pi */
  y = AVL_mul_ps(x, *(AVL_FloatVectorType*)_ps512_cephes_FOPI);

  /* store the integer part of y in imm2 */
  imm2 = AVL_cvttps_epi32(y);

  /* j=(j+1) & (~1) (see the cephes sources) */
  imm2 = AVL_add_epi32(imm2, *(AVL_IntegerVectorType*)_pi32_512_1);
  imm2 = AVL_and_si(imm2, *(AVL_IntegerVectorType*)_pi32_512_inv1);

  y = AVL_cvtepi32_ps(imm2);
  imm4 = imm2;

  /* get the swap sign flag for the sine */
  imm0 = AVL_and_si(imm2, *(AVL_IntegerVectorType*)_pi32_512_4);
  imm0 = AVL_slli_epi32(imm0, 29);
  //AVL_FloatVectorType swap_sign_bit_sin = AVL_castsi_ps(imm0);

  /* get the polynom selection mask for the sine*/
  imm2 = AVL_and_si(imm2, *(AVL_IntegerVectorType*)_pi32_512_2);
  AVL_MaskType imm22 = AVL_cmpeq_epi32_mask(imm2, *(AVL_IntegerVectorType*)_pi32_512_0);
  imm2 = AVL_mask_blend_epi32(imm22, *(AVL_IntegerVectorType*)_pi32_512_0, *(AVL_IntegerVectorType*)_pi32_512_0xffffffff);

  //AVL_FloatVectorType poly_mask = AVL_castsi_ps(imm2);

  AVL_FloatVectorType swap_sign_bit_sin = AVL_castsi_ps(imm0);
  AVL_FloatVectorType poly_mask = AVL_castsi_ps(imm2);

  /* The magic pass: "Extended precision modular arithmetic"
     x = ((x - y * DP1) - y * DP2) - y * DP3; */
  xmm1 = *(AVL_FloatVectorType*)_ps512_minus_cephes_DP1;
  xmm2 = *(AVL_FloatVectorType*)_ps512_minus_cephes_DP2;
  xmm3 = *(AVL_FloatVectorType*)_ps512_minus_cephes_DP3;
  xmm1 = AVL_mul_ps(y, xmm1);
  xmm2 = AVL_mul_ps(y, xmm2);
  xmm3 = AVL_mul_ps(y, xmm3);
  x = AVL_add_ps(x, xmm1);
  x = AVL_add_ps(x, xmm2);
  x = AVL_add_ps(x, xmm3);

  imm4 = AVL_sub_epi32(imm4, *(AVL_IntegerVectorType*)_pi32_512_2);
  imm4 = AVL_andnot_si(imm4, *(AVL_IntegerVectorType*)_pi32_512_4);
  imm4 = AVL_slli_epi32(imm4, 29);

  AVL_FloatVectorType sign_bit_cos = AVL_castsi_ps(imm4);

  sign_bit_sin = AVL_castsi_ps(AVL_xor_si(AVL_castps_si(sign_bit_sin), AVL_castps_si(swap_sign_bit_sin)));

  /* Evaluate the first polynom  (0 <= x <= Pi/4) */
  AVL_FloatVectorType z = AVL_mul_ps(x, x);
  y = *(AVL_FloatVectorType*)_ps512_coscof_p0;

  y = AVL_mul_ps(y, z);
  y = AVL_add_ps(y, *(AVL_FloatVectorType*)_ps512_coscof_p1);
  y = AVL_mul_ps(y, z);
  y = AVL_add_ps(y, *(AVL_FloatVectorType*)_ps512_coscof_p2);
  y = AVL_mul_ps(y, z);
  y = AVL_mul_ps(y, z);
  AVL_FloatVectorType tmp = AVL_mul_ps(z, *(AVL_FloatVectorType*)_ps512_0p5);
  y = AVL_sub_ps(y, tmp);
  y = AVL_add_ps(y, *(AVL_FloatVectorType*)_ps512_1);

  /* Evaluate the second polynom  (Pi/4 <= x <= 0) */

  AVL_FloatVectorType y2 = *(AVL_FloatVectorType*)_ps512_sincof_p0;
  y2 = AVL_mul_ps(y2, z);
  y2 = AVL_add_ps(y2, *(AVL_FloatVectorType*)_ps512_sincof_p1);
  y2 = AVL_mul_ps(y2, z);
  y2 = AVL_add_ps(y2, *(AVL_FloatVectorType*)_ps512_sincof_p2);
  y2 = AVL_mul_ps(y2, z);
  y2 = AVL_mul_ps(y2, x);
  y2 = AVL_add_ps(y2, x);

  /* select the correct result from the two polynoms */
  xmm3 = poly_mask;
  AVL_FloatVectorType ysin2 = AVL_castsi_ps(AVL_and_si(AVL_castps_si(xmm3), AVL_castps_si(y2)));
  AVL_FloatVectorType ysin1 = AVL_castsi_ps(AVL_andnot_si(AVL_castps_si(xmm3), AVL_castps_si(y)));
  y2 = AVL_sub_ps(y2, ysin2);
  y = AVL_sub_ps(y, ysin1);

  xmm1 = AVL_add_ps(ysin1, ysin2);
  xmm2 = AVL_add_ps(y, y2);

  /* update the sign */
  *s = AVL_castsi_ps(AVL_xor_si(AVL_castps_si(xmm1), AVL_castps_si(sign_bit_sin)));
  *c = AVL_castsi_ps(AVL_xor_si(AVL_castps_si(xmm2), AVL_castps_si(sign_bit_cos)));
}


int roundUp(std::size_t numToRound, std::size_t multiple)
{
  if (multiple == 0)
    return numToRound;

  int remainder = numToRound % multiple;
  if (remainder == 0)
    return numToRound;

  return numToRound + multiple - remainder;
}


AVXVectorMask::AVXVectorMask(const std::size_t& size): _size(size) {
  if (size == 0) {
    throw std::runtime_error("Creating AVXVectorMask with size 0 is not allowed");
  }

  _data = reinterpret_cast<AVL_MaskType*>(_mm_malloc(roundUp(_size / AVL_VectorSize, AVL_VectorSize) * sizeof(AVL_MaskType), AVL_Alignment));

}
AVXVectorMask::~AVXVectorMask() {
  _mm_free(_data);
}
void AVXVectorMask::setData(const AVXVectorMask* data) {
  if (_size != data->_size ) {
    throw std::runtime_error("Not the same size in setData");
  }
  memcpy(_data, data->_data, roundUp(_size / AVL_VectorSize, AVL_VectorSize) * sizeof(AVL_MaskType));

}
void AVXVectorMask::setData(const bool& data) {
  if (data) {
    memset ( _data, 0xff, roundUp(_size / AVL_VectorSize, AVL_VectorSize) * sizeof(AVL_MaskType) );
  }
  else {
    memset ( _data, 0, roundUp(_size / AVL_VectorSize, AVL_VectorSize) * sizeof(AVL_MaskType) );
  }

}


void AVXVectorMask::OR(const AVXVectorMask* result, const AVXVectorMask* a, const AVXVectorMask* b) {
  if (result->_size != a->_size || result->_size != b->_size) {
    throw std::runtime_error("Not the same size in OR");
  }
  std::size_t iterations = result->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_MaskType a_data = AVL_load_mask(a->_data + i);
    AVL_MaskType b_data = AVL_load_mask(b->_data + i);
    AVL_MaskType result_data = AVL_kor(a_data, b_data);
    AVL_store_mask(result->_data + i, result_data);
  }

}

void AVXVectorMask::AND(const AVXVectorMask* result, const AVXVectorMask* a, const AVXVectorMask* b) {
  if (result->_size != a->_size || result->_size != b->_size) {
    throw std::runtime_error("Not the same size in OR");
  }
  std::size_t iterations = result->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_MaskType a_data = AVL_load_mask(a->_data + i);
    AVL_MaskType b_data = AVL_load_mask(b->_data + i);
    AVL_MaskType result_data = AVL_kand(a_data, b_data);
    AVL_store_mask(result->_data + i, result_data);
  }

}

void AVXVectorMask::NOT(const AVXVectorMask* result, const AVXVectorMask* a) {
  if (result->_size != a->_size ) {
    throw std::runtime_error("Not the same size in NOT");
  }
  std::size_t iterations = result->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_MaskType a_data = AVL_load_mask(a->_data + i);
    AVL_MaskType result_data = AVL_knot(a_data);
    AVL_store_mask(result->_data + i, result_data);
  }

}


std::size_t AVXVectorMask::getIterations() const {
  std::size_t it = _size / AVL_VectorSize + 1;
  if (_size % AVL_VectorSize == 0) {
    it -= 1;
  }
  return it;

}
AVL_MaskType AVXVectorF::getLastMask() const{
  uint8_t diff=getIterations()*AVL_VectorSize-_size;
  AVL_MaskType mask=AVL_cvtu32_mask(0);
  mask=AVL_knot(mask);
  //AVL_MaskType mask=0xffff;
  mask=mask>>diff;
  return mask;

}

AVXVectorF::AVXVectorF(const std::size_t& size): _size(size) {
  if (size == 0) {
    throw std::runtime_error("Creating AVXVector with size 0 is not allowed");
  }

  _data = reinterpret_cast<float*>(_mm_malloc(roundUp(_size, AVL_VectorSize) * sizeof(float), AVL_Alignment));

}
AVXVectorF::~AVXVectorF() {
  _mm_free(_data);
}

void AVXVectorF::setData(const float& data) {
  std::size_t iterations = getIterations();
  AVL_FloatVectorType constant = AVL_set1_ps(data);
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_store_ps(_data + i * AVL_VectorSize, constant);
  }

}

void AVXVectorF::setData(const float* data) {
  memcpy(_data, data, _size * sizeof(float));
}

void AVXVectorF::setData(const int16_t* data) {
  std::size_t iterations = getIterations();
  AVL_MaskType mask=AVL_cvtu32_mask(0);
  mask=AVL_knot(mask);
//  AVL_MaskType mask = 0xffff;
  for (std::size_t i = 0; i < iterations; i++) {
    //load 16 signed 16bit integers
    //__m256i t1= _mm256_loadu_epi16 (data+i*16);
    AVL_IntegerVectorConversionType t1 = AVL_maskz_loadu_epi16 (mask, data + i * 16);
    //convert to 32bit signed integers
    AVL_IntegerVectorType t2 = AVL_cvtepi16_epi32 (t1);
    //convert to float
    AVL_FloatVectorType t3  = AVL_cvtepi32_ps  (t2);
    //store
    AVL_store_ps(_data + i * AVL_VectorSize, t3);
  }


}
void AVXVectorF::setData(const std::vector<float>& data) {
  if(data.size()!=_size){
    throw std::runtime_error("Not the same size in setData");
  }
  std::copy(data.begin(), data.end(), _data);
}

void AVXVectorF::setData(const AVXVectorF* data) {
  if (_size != data->_size ) {
    throw std::runtime_error("Not the same size in setData");
  }
  memcpy(_data, data->_data, _size * (sizeof(float)));
}
void AVXVectorF::setData(const AVXVectorMask* mask, const float& a) {
  if (_size != mask->_size ) {
    throw std::runtime_error("Not the same size in setData");
  }
  std::size_t iterations = getIterations();
  AVL_FloatVectorType constant = AVL_set1_ps(a);
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_MaskType mask_data = AVL_load_mask(mask->_data + i);
    AVL_mask_store_ps(_data + i * AVL_VectorSize, mask_data, constant);

  }

}

void AVXVectorF::print() {
  for (std::size_t i = 0; i < _size; i++) {
    std::cout << std::setprecision(10) << *(_data + i) << ",";
  }
  std::cout << std::endl;
}

std::size_t AVXVectorF::getIterations() const {
  std::size_t it = _size / AVL_VectorSize + 1;
  if (_size % AVL_VectorSize == 0) {
    it -= 1;
  }
  return it;

}

void AVXVectorF::add(const AVXVectorF* result, const AVXVectorF* a, const AVXVectorF* b) {
  if (result->_size != a->_size || result->_size != b->_size) {
    throw std::runtime_error("Not the same size in add");
  }
  std::size_t iterations = result->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType b_data = AVL_load_ps (b->_data + i * AVL_VectorSize);
    AVL_FloatVectorType result_data = AVL_add_ps(a_data, b_data);
    AVL_store_ps(result->_data + i * AVL_VectorSize, result_data);

  }


}

void AVXVectorF::add(const AVXVectorF* result, const AVXVectorF* a, const float& b) {
  if (result->_size != a->_size ) {
    throw std::runtime_error("Not the same size in add");
  }
  std::size_t iterations = result->getIterations();
  AVL_FloatVectorType constant = AVL_set1_ps(b);
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType result_data = AVL_add_ps(a_data, constant);
    AVL_store_ps(result->_data + i * AVL_VectorSize, result_data);

  }
}

void AVXVectorF::add_mask(const AVXVectorF* result, const AVXVectorMask* mask, const AVXVectorF* a, const float& b) {
  if (result->_size != a->_size ) {
    throw std::runtime_error("Not the same size in add");
  }
  std::size_t iterations = result->getIterations();
  AVL_FloatVectorType constant = AVL_set1_ps(b);
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_MaskType mask_data = AVL_load_mask(mask->_data + i);
    AVL_FloatVectorType result_data = AVL_load_ps (result->_data + i * AVL_VectorSize);
    result_data = AVL_mask_add_ps(result_data, mask_data, a_data, constant);
    AVL_store_ps(result->_data + i * AVL_VectorSize, result_data);

  }
}



void AVXVectorF::sub(const AVXVectorF* result, const AVXVectorF* a, const AVXVectorF* b) {
  if (result->_size != a->_size || result->_size != b->_size) {
    throw std::runtime_error("Not the same size in add");
  }
  std::size_t iterations = result->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType b_data = AVL_load_ps (b->_data + i * AVL_VectorSize);
    AVL_FloatVectorType result_data = AVL_sub_ps(a_data, b_data);
    AVL_store_ps(result->_data + i * AVL_VectorSize, result_data);

  }
}

void AVXVectorF::mult(const AVXVectorF* result, const AVXVectorF* a, const AVXVectorF* b) {
  if (result->_size != a->_size || result->_size != b->_size) {
    throw std::runtime_error("Not the same size in mult");
  }
  std::size_t iterations = result->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType b_data = AVL_load_ps (b->_data + i * AVL_VectorSize);
    AVL_FloatVectorType result_data = AVL_mul_ps(a_data, b_data);
    AVL_store_ps(result->_data + i * AVL_VectorSize, result_data);

  }
}

void AVXVectorF::mult(const AVXVectorF* result, const AVXVectorF* a, const float& b) {
  if (result->_size != a->_size ) {
    throw std::runtime_error("Not the same size in add");
  }
  std::size_t iterations = result->getIterations();
  AVL_FloatVectorType constant = AVL_set1_ps(b);
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType result_data = AVL_mul_ps(a_data, constant);
    AVL_store_ps(result->_data + i * AVL_VectorSize, result_data);

  }
}

void AVXVectorF::div(const AVXVectorF* result, const AVXVectorF* a, const AVXVectorF* b) {
  if (result->_size != a->_size || result->_size != b->_size) {
    throw std::runtime_error("Not the same size in mult");
  }
  std::size_t iterations = result->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType b_data = AVL_load_ps (b->_data + i * AVL_VectorSize);
    AVL_FloatVectorType result_data = AVL_div_ps(a_data, b_data);

    AVL_store_ps(result->_data + i * AVL_VectorSize, result_data);

  }
}

void AVXVectorF::divWithInf(const AVXVectorF* result, const AVXVectorF* a, const AVXVectorF* b) {
  if (result->_size != a->_size || result->_size != b->_size) {
    throw std::runtime_error("Not the same size in mult");
  }
  AVL_FloatVectorType zero = AVL_set1_ps(0.0f);
  const AVL_FloatVectorType INF = AVL_set1_ps(std::numeric_limits<float>::infinity());
  const AVL_FloatVectorType NEG_INF = AVL_set1_ps(-std::numeric_limits<float>::infinity());
  std::size_t iterations = result->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType b_data = AVL_load_ps (b->_data + i * AVL_VectorSize);
    AVL_FloatVectorType result_data = AVL_div_ps(a_data, b_data);

    AVL_MaskType div_with_zero = AVL_cmp_ps_mask(b_data, zero, _CMP_EQ_OQ);
    AVL_MaskType greater_than_zero = AVL_cmp_ps_mask(a_data, zero, _CMP_GT_OQ);
    AVL_MaskType less_than_zero = AVL_cmp_ps_mask(a_data, zero, _CMP_LT_OQ);

    AVL_MaskType inf = AVL_kand(div_with_zero, greater_than_zero);
    AVL_MaskType neg_inf = AVL_kand(div_with_zero, less_than_zero);

    result_data = AVL_mask_mov_ps(result_data, inf, INF);
    result_data = AVL_mask_mov_ps(result_data, neg_inf, NEG_INF);

    AVL_store_ps(result->_data + i * AVL_VectorSize, result_data);

  }
}

void AVXVectorF::div(const AVXVectorF* result, const AVXVectorF* a, const float& b) {
  if (result->_size != a->_size ) {
    throw std::runtime_error("Not the same size in add");
  }
  std::size_t iterations = result->getIterations();
  AVL_FloatVectorType constant = AVL_set1_ps(b);
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType result_data = AVL_div_ps(a_data, constant);
    AVL_store_ps(result->_data + i * AVL_VectorSize, result_data);

  }
}



void AVXVectorF::fmadd(const AVXVectorF* result, const AVXVectorF* a, const AVXVectorF* b, const AVXVectorF* c) {
  if (result->_size != a->_size || result->_size != b->_size || result->_size != c->_size) {
    throw std::runtime_error("Not the same size in fmadd");
  }
  std::size_t iterations = result->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType b_data = AVL_load_ps (b->_data + i * AVL_VectorSize);
    AVL_FloatVectorType c_data = AVL_load_ps (c->_data + i * AVL_VectorSize);
    AVL_FloatVectorType result_data = AVL_fmadd_ps(a_data, b_data, c_data);
    AVL_store_ps(result->_data + i * AVL_VectorSize, result_data);
  }
}

void AVXVectorF::fmadd(const AVXVectorF* result, const AVXVectorF* a, const float& b, const AVXVectorF* c) {
  if (result->_size != a->_size ||  result->_size != c->_size) {
    throw std::runtime_error("Not the same size in fmadd");
  }
  std::size_t iterations = result->getIterations();
  AVL_FloatVectorType constant = AVL_set1_ps(b);
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType c_data = AVL_load_ps (c->_data + i * AVL_VectorSize);
    AVL_FloatVectorType result_data = AVL_fmadd_ps(a_data, constant, c_data);
    AVL_store_ps(result->_data + i * AVL_VectorSize, result_data);

  }
}

void AVXVectorF::fmsub(const AVXVectorF* result, const AVXVectorF* a, const AVXVectorF* b, const AVXVectorF* c) {
  if (result->_size != a->_size || result->_size != b->_size || result->_size != c->_size) {
    throw std::runtime_error("Not the same size in fmadd");
  }
  std::size_t iterations = result->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType b_data = AVL_load_ps (b->_data + i * AVL_VectorSize);
    AVL_FloatVectorType c_data = AVL_load_ps (c->_data + i * AVL_VectorSize);
    AVL_FloatVectorType result_data = AVL_fmsub_ps(a_data, b_data, c_data);
    AVL_store_ps(result->_data + i * AVL_VectorSize, result_data);

  }
}



void AVXVectorF::ln(const AVXVectorF* result, const AVXVectorF* a) {
  if (result->_size != a->_size ) {
    throw std::runtime_error("Not the same size in ln");
  }
  std::size_t iterations = result->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType result_data = log512_ps(a_data);
    AVL_store_ps(result->_data + i * AVL_VectorSize, result_data);

  }
}

void AVXVectorF::exp(const AVXVectorF* result, const AVXVectorF* a) {
  if (result->_size != a->_size ) {
    throw std::runtime_error("Not the same size in ln");
  }
  std::size_t iterations = result->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType result_data = exp512_ps(a_data);
    AVL_store_ps(result->_data + i * AVL_VectorSize, result_data);

  }
}

void AVXVectorF::sqrt(const AVXVectorF* result, const AVXVectorF* a) {
  if (result->_size != a->_size ) {
    throw std::runtime_error("Not the same size in ln");
  }
  std::size_t iterations = result->getIterations();
  AVL_MaskType temp_mask=AVL_cvtu32_mask(0);
  temp_mask=AVL_knot(temp_mask);
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType result_data = AVL_maskz_rsqrt14_ps(temp_mask,a_data);
    result_data = AVL_rcp14_ps(result_data);
    AVL_store_ps(result->_data + i * AVL_VectorSize, result_data);

  }
}

void AVXVectorF::complexDiv(const AVXVectorF* resultReal, const AVXVectorF* resultImag, const AVXVectorF* a, const AVXVectorF* b, const AVXVectorF* c, const AVXVectorF* d) {
  if (resultReal->_size != resultImag->_size || resultReal->_size != a->_size || resultReal->_size != b->_size || resultReal->_size != c->_size || resultReal->_size != d->_size ) {
    throw std::runtime_error("Not the same size in complexDiv");
  }
  AVL_FloatVectorType zero = AVL_set1_ps(0.0f);
  const AVL_FloatVectorType INF = AVL_set1_ps(std::numeric_limits<float>::infinity());
  const AVL_FloatVectorType NEG_INF = AVL_set1_ps(-std::numeric_limits<float>::infinity());
  std::size_t iterations = resultReal->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType b_data = AVL_load_ps (b->_data + i * AVL_VectorSize);
    AVL_FloatVectorType c_data = AVL_load_ps (c->_data + i * AVL_VectorSize);
    AVL_FloatVectorType d_data = AVL_load_ps (d->_data + i * AVL_VectorSize);

    AVL_FloatVectorType tempReal = AVL_mul_ps(b_data, d_data);
    tempReal = AVL_fmadd_ps(a_data, c_data, tempReal);

    AVL_FloatVectorType tempImag = AVL_mul_ps(a_data, d_data);
    tempImag = AVL_fmsub_ps(b_data, c_data, tempImag);

    AVL_FloatVectorType tempDiv = AVL_mul_ps(d_data, d_data);
    tempDiv = AVL_fmadd_ps(c_data, c_data, tempDiv);
    AVL_FloatVectorType real = AVL_div_ps(tempReal, tempDiv);
    AVL_FloatVectorType imag = AVL_div_ps(tempImag, tempDiv);

    AVL_MaskType div_with_zero = AVL_cmp_ps_mask(tempDiv, zero, _CMP_EQ_OQ);
   // std::cout<<"div with zero:"<<std::endl;
    //print(div_with_zero);
    AVL_MaskType greater_than_zero = AVL_cmp_ps_mask(a_data, zero, _CMP_GT_OQ);
    //std::cout<<"real greater than zero:"<<std::endl;
    //print(greater_than_zero);
    AVL_MaskType less_than_zero = AVL_cmp_ps_mask(a_data, zero, _CMP_LT_OQ);
    //std::cout<<"real less than zero:"<<std::endl;
    //print(less_than_zero);
    AVL_MaskType inf = AVL_kand(div_with_zero, greater_than_zero);
    AVL_MaskType neg_inf = AVL_kand(div_with_zero, less_than_zero);

    real = AVL_mask_mov_ps(real, inf, INF);
    real = AVL_mask_mov_ps(real, neg_inf, NEG_INF);

     greater_than_zero = AVL_cmp_ps_mask(b_data, zero, _CMP_GT_OQ);
     less_than_zero = AVL_cmp_ps_mask(b_data, zero, _CMP_LT_OQ);

     inf = AVL_kand(div_with_zero, greater_than_zero);
     neg_inf = AVL_kand(div_with_zero, less_than_zero);

    imag = AVL_mask_mov_ps(imag, inf, INF);
    imag = AVL_mask_mov_ps(imag, neg_inf, NEG_INF);

    AVL_store_ps(resultReal->_data + i * AVL_VectorSize, real);
    AVL_store_ps(resultImag->_data + i * AVL_VectorSize, imag);

  }
}



void AVXVectorF::atan2(const AVXVectorF* result, const AVXVectorF* a, const AVXVectorF* b) {
  if (result->_size != a->_size || result->_size != b->_size) {
    throw std::runtime_error("Not the same size in mult");
  }
  std::size_t iterations = result->getIterations();
  AVL_FloatVectorType ps_0 = AVL_set1_ps(0.0f);
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType y = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType x = AVL_load_ps (b->_data + i * AVL_VectorSize);
    AVL_FloatVectorType _x = AVL_div_ps(y, x);
  // std::cout<<"----------own-----------"<<std::endl;
  //  std::cout << "y/x:" << std::endl;
  //  print(_x);
    //make argument postive and save the sign
    AVL_MaskType sign_mask = AVL_cmp_ps_mask(_x, ps_0, _CMP_LE_OQ);
    _x = AVL_abs_ps(_x);

    //mask = _x>2.414.....
    AVL_MaskType range_reduction0 = AVL_cmp_ps_mask(_x, AVL_set1_ps(2.414213562373095), _CMP_GT_OQ);
    //mask = _x>0.414....
    AVL_MaskType range_reduction1 = AVL_cmp_ps_mask(_x, AVL_set1_ps(0.4142135623730950), _CMP_GT_OQ);
    //mask =_x>0.414 and _x<=2.414
    AVL_MaskType range_reduction2 = AVL_knot(range_reduction0);
    range_reduction2 = AVL_kand(range_reduction1, range_reduction2);

    //x_inv_neg=-1/_x
    AVL_FloatVectorType x_inv_neg = AVL_rcp14_ps (_x);
    x_inv_neg = AVL_sub_ps(AVL_set1_ps(0.0f), x_inv_neg);

    //if _x>2.414... then _x=-1/_x
    _x = AVL_mask_mov_ps (_x, range_reduction0, x_inv_neg);
    AVL_FloatVectorType _y = AVL_set1_ps(1.5707963267948966192);


    AVL_FloatVectorType x_minus_one = _x;
    x_minus_one = AVL_sub_ps(x_minus_one, AVL_set1_ps(1.0f));
    AVL_FloatVectorType x_plus_one = _x;
    x_plus_one = AVL_add_ps(x_plus_one, AVL_set1_ps(1.0f));
    //if _x>0.414 and _x<=2.414 then _x= ( _x - 1.0 ) / ( _x + 1.0 )
    _x = AVL_mask_div_ps (_x, range_reduction2, x_minus_one, x_plus_one);
    //if _x>0.414 and _x<=2.414 then _y=0.7853.....
    _y = AVL_mask_mov_ps(_y, range_reduction2, AVL_set1_ps(0.7853981633974483096));
    range_reduction1 = AVL_knot(range_reduction1);
    //if _x<=0.414.... _y=0
    _y = AVL_mask_mov_ps(_y, range_reduction1, AVL_set1_ps(0.0f));
   // std::cout << "x after range reduction:" << std::endl;
   // print(_x);
   // std::cout << "y after range reduction:" << std::endl;
    //print(_y);


    AVL_FloatVectorType z = AVL_mul_ps(_x, _x);
    //std::cout << "z:" << std::endl;
    //print(z);
    AVL_FloatVectorType temp = AVL_fmsub_ps(AVL_set1_ps( 8.05374449538e-2), z, AVL_set1_ps( 1.38776856032E-1));
    temp = AVL_fmadd_ps(temp, z, AVL_set1_ps( 1.99777106478E-1));
    temp = AVL_fmsub_ps(temp, z, AVL_set1_ps( 3.33329491539E-1));
    temp = AVL_mul_ps(temp, z);
    temp = AVL_fmadd_ps(temp, _x, _x);
    _y = AVL_add_ps(_y, temp);
    //std::cout << "y after approximation:" << std::endl;
    //print(_y);
    _y = AVL_mask_sub_ps (_y, sign_mask, AVL_set1_ps(0.0f), _y);


    //y>0 x=0
    AVL_MaskType temp_mask = AVL_cmp_ps_mask(y, ps_0, _CMP_GT_OQ);
    AVL_MaskType mode1 = AVL_cmp_ps_mask(x, ps_0, _CMP_EQ_OQ);
    mode1 = AVL_kand(temp_mask, mode1);

    //y>0 x<0
    temp_mask = AVL_cmp_ps_mask(y, ps_0, _CMP_GT_OQ);
    AVL_MaskType mode2 = AVL_cmp_ps_mask(x, ps_0, _CMP_LT_OQ);
    mode2 = AVL_kand(temp_mask, mode2);

    //y=0 x>0
    temp_mask = AVL_cmp_ps_mask(y, ps_0, _CMP_EQ_OQ);
    AVL_MaskType mode3 = AVL_cmp_ps_mask(x, ps_0, _CMP_GT_OQ);
    mode3 = AVL_kand(temp_mask, mode3);

    //y=0 x=0
    temp_mask = AVL_cmp_ps_mask(y, ps_0, _CMP_EQ_OQ);
    AVL_MaskType mode4 = AVL_cmp_ps_mask(x, ps_0, _CMP_EQ_OQ);
    mode4 = AVL_kand(temp_mask, mode4);

    //y=0 x<0
    temp_mask = AVL_cmp_ps_mask(y, ps_0, _CMP_EQ_OQ);
    AVL_MaskType mode5 = AVL_cmp_ps_mask(x, ps_0, _CMP_LT_OQ);
    mode5 = AVL_kand(temp_mask, mode5);

    //y<0 x=0
    temp_mask = AVL_cmp_ps_mask(y, ps_0, _CMP_LT_OQ);
    AVL_MaskType mode7 = AVL_cmp_ps_mask(x, ps_0, _CMP_EQ_OQ);
    mode7 = AVL_kand(temp_mask, mode7);

    //y<0 x<0
    temp_mask = AVL_cmp_ps_mask(y, ps_0, _CMP_LT_OQ);
    AVL_MaskType mode8 = AVL_cmp_ps_mask(x, ps_0, _CMP_LT_OQ);
    mode8 = AVL_kand(temp_mask, mode8);

    AVL_FloatVectorType w = AVL_set1_ps(0.0f);
    w = AVL_mask_mov_ps(w, mode2, AVL_set1_ps(3.141592653589793238));
    w = AVL_mask_mov_ps(w, mode8, AVL_set1_ps(-3.141592653589793238));
   // std::cout << "w:" << std::endl;
   // print(w);
    _y = AVL_mask_mov_ps(_y, mode1, AVL_set1_ps(1.5707963267948966192));
    _y = AVL_mask_mov_ps(_y, mode3, AVL_set1_ps(0.0f));
    _y = AVL_mask_mov_ps(_y, mode4, AVL_set1_ps(0.0f));
    _y = AVL_mask_mov_ps(_y, mode5, AVL_set1_ps(3.141592653589793238));
    _y = AVL_mask_mov_ps(_y, mode7, AVL_set1_ps(-1.5707963267948966192));
   // std::cout << "x:" << std::endl;
   // print(_y);
    _y = AVL_add_ps(_y, w);
    AVL_store_ps(result->_data + i * AVL_VectorSize, _y);

  }
}

std::size_t AVXVectorF::count(const AVXVectorF * a,const float& b){
  AVL_FloatVectorType constant = AVL_set1_ps(b);
  std::size_t temp=0;
  std::size_t iterations = a->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_MaskType result_data = AVL_cmp_ps_mask(a_data, constant, _CMP_EQ_OQ);
    for(std::size_t i=0;i<AVL_VectorSize;i++){
      if(result_data>>i & 0x1){
        temp++;
      }
    }
  }
  return temp;

}



float AVXVectorF::sum(const AVXVectorF* a) {
  AVL_FloatVectorType result = AVL_set1_ps(0.0f);
  std::size_t iterations = a->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    if(i!=iterations-1){
      result = AVL_add_ps(a_data, result);
    }
    else{
      result= AVL_mask_add_ps(result,a->getLastMask(),result,a_data);
    }
  }
  float temp_buffer[AVL_VectorSize];
  float to_return =0;
  AVL_store_ps(temp_buffer,result);
  for(std::size_t i =0;i<AVL_VectorSize;i++){
    to_return+=temp_buffer[i];
  }
  return to_return;
  //return _mm512_reduce_add_ps(result);

}
  float AVXVectorF::min(const AVXVectorF * a){
  AVL_FloatVectorType result = AVL_load_ps (a->_data);
  std::size_t iterations = a->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    if(i!=iterations-1){
      result = AVL_min_ps(a_data, result);
    }
    else{
      result= AVL_mask_min_ps(result,a->getLastMask(),result,a_data);
    }
  }
  float temp_buffer[AVL_VectorSize];
  AVL_store_ps(temp_buffer,result);
  float to_return =temp_buffer[0];
  for(std::size_t i =1;i<AVL_VectorSize;i++){
    to_return=std::min(temp_buffer[i],to_return);
  }
  return to_return;
  //return _mm512_reduce_min_ps(result);
  }

  float AVXVectorF::minNotZero(const AVXVectorF * a){
  AVL_FloatVectorType result = AVL_set1_ps(std::numeric_limits<float>::infinity());
  AVL_FloatVectorType constant = AVL_set1_ps (0.0f);
  std::size_t iterations = a->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {

    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_MaskType mask = AVL_cmp_ps_mask(a_data, constant, _CMP_NEQ_OQ);
    if(i!=iterations-1){
      result = AVL_mask_min_ps(result,mask,a_data, result);
    }
    else{
      mask=AVL_kand(mask, a->getLastMask());
      result= AVL_mask_min_ps(result,mask,result,a_data);
    }
  }
  float temp_buffer[AVL_VectorSize];
  AVL_store_ps(temp_buffer,result);
  float to_return =temp_buffer[0];
  for(std::size_t i =1;i<AVL_VectorSize;i++){
    to_return=std::min(temp_buffer[i],to_return);
  }
  return to_return;
  //return _mm512_reduce_min_ps(result);

  }

  float AVXVectorF::max(const AVXVectorF * a){
  AVL_FloatVectorType result = AVL_load_ps (a->_data);
  std::size_t iterations = a->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    if(i!=iterations-1){
      result = AVL_max_ps(a_data, result);
    }
    else{
      result= AVL_mask_max_ps(result,a->getLastMask(),result,a_data);
    }
  }
  float temp_buffer[AVL_VectorSize];
  AVL_store_ps(temp_buffer,result);
  float to_return =temp_buffer[0];
  for(std::size_t i =1;i<AVL_VectorSize;i++){
    to_return=std::max(temp_buffer[i],to_return);
  }
  return to_return;
//  return _mm512_reduce_max_ps(result);
  }

void AVXVectorF::rmInfNan(const AVXVectorF* a) {
  const AVL_FloatVectorType SIGN_MASK = AVL_set1_ps(-0.0);
  const AVL_FloatVectorType INF = AVL_set1_ps(std::numeric_limits<float>::infinity());
  std::size_t iterations = a->getIterations();
  AVL_FloatVectorType zeros = AVL_setzero_ps ();
  for (std::size_t i = 0; i < iterations; i++) {
    //load data
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    //mask out sign bit
    AVL_FloatVectorType sign_mask = AVL_andnot_ps(SIGN_MASK, a_data);
    //create mask for infinity
    AVL_MaskType mask = AVL_cmp_ps_mask(sign_mask, INF, _CMP_EQ_OQ);
    //create mask for nan
    AVL_MaskType nan_mask = AVL_cmp_ps_mask (a_data, a_data, _CMP_NEQ_UQ);
    //or both masks
    mask =   AVL_kor (nan_mask, mask);
    //replace infinity with zeros
    a_data = AVL_mask_blend_ps(mask, a_data, zeros);
    //store result
    AVL_store_ps(a->_data + i * AVL_VectorSize, a_data);
  }
}

void AVXVectorF::min(const AVXVectorF* result, const AVXVectorF* a, const AVXVectorF* b) {
  if (result->_size != a->_size || result->_size != b->_size) {
    throw std::runtime_error("Not the same size in min");
  }
  std::size_t iterations = result->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType b_data = AVL_load_ps (b->_data + i * AVL_VectorSize);
    AVL_FloatVectorType result_data = AVL_min_ps(a_data, b_data);
    AVL_store_ps(result->_data + i * AVL_VectorSize, result_data);

  }


}

void AVXVectorF::max(const AVXVectorF* result, const AVXVectorF* a, const AVXVectorF* b) {
  if (result->_size != a->_size || result->_size != b->_size) {
    throw std::runtime_error("Not the same size in max");
  }
  std::size_t iterations = result->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType b_data = AVL_load_ps (b->_data + i * AVL_VectorSize);
    AVL_FloatVectorType result_data = AVL_max_ps(a_data, b_data);
    AVL_store_ps(result->_data + i * AVL_VectorSize, result_data);

  }


}

void AVXVectorF::equal(const AVXVectorMask* result, const AVXVectorF* a, const AVXVectorF* b) {
  if (result->_size != a->_size || result->_size != b->_size) {
    throw std::runtime_error("Not the same size in max");
  }
  std::size_t iterations = a->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType b_data = AVL_load_ps (b->_data + i * AVL_VectorSize);
    AVL_MaskType result_data = AVL_cmp_ps_mask(a_data, b_data, _CMP_EQ_OQ);
    AVL_store_mask(result->_data + i, result_data);

  }
}

void AVXVectorF::equal(const AVXVectorMask* result, const AVXVectorF* a, const float& b) {
  if (result->_size != a->_size) {
    throw std::runtime_error("Not the same size in max");
  }
  AVL_FloatVectorType constant = AVL_set1_ps(b);
  std::size_t iterations = a->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_MaskType result_data = AVL_cmp_ps_mask(a_data, constant, _CMP_EQ_OQ);
    AVL_store_mask(result->_data + i, result_data);

  }
}

void AVXVectorF::greater_than(const AVXVectorMask* result, const AVXVectorF* a, const AVXVectorF* b) {
  if (result->_size != a->_size || result->_size != b->_size) {
    throw std::runtime_error("Not the same size in max");
  }
  std::size_t iterations = a->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType b_data = AVL_load_ps (b->_data + i * AVL_VectorSize);
    AVL_MaskType result_data = AVL_cmp_ps_mask(a_data, b_data, _CMP_GT_OQ);
    AVL_store_mask(result->_data + i, result_data);

  }
}

void AVXVectorF::greater_than(const AVXVectorMask* result, const AVXVectorF* a, const float& b) {
  if (result->_size != a->_size ) {
    throw std::runtime_error("Not the same size in max");
  }
  AVL_FloatVectorType constant = AVL_set1_ps(b);
  std::size_t iterations = a->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_MaskType result_data = AVL_cmp_ps_mask(a_data, constant, _CMP_GT_OQ);
    AVL_store_mask(result->_data + i, result_data);

  }
}

void AVXVectorF::greater_equal_than(const AVXVectorMask* result, const AVXVectorF* a, const AVXVectorF* b) {
  if (result->_size != a->_size || result->_size != b->_size) {
    throw std::runtime_error("Not the same size in max");
  }
  std::size_t iterations = a->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType b_data = AVL_load_ps (b->_data + i * AVL_VectorSize);
    AVL_MaskType result_data = AVL_cmp_ps_mask(a_data, b_data, _CMP_GE_OQ);
    AVL_store_mask(result->_data + i, result_data);

  }
}

void AVXVectorF::less_than(const AVXVectorMask* result, const AVXVectorF* a, const AVXVectorF* b) {
  if (result->_size != a->_size || result->_size != b->_size) {
    throw std::runtime_error("Not the same size in max");
  }
  std::size_t iterations = a->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType b_data = AVL_load_ps (b->_data + i * AVL_VectorSize);
    AVL_MaskType result_data = AVL_cmp_ps_mask(a_data, b_data, _CMP_LT_OQ);
    AVL_store_mask(result->_data + i, result_data);

  }
}

void AVXVectorF::less_equal_than(const AVXVectorMask* result, const AVXVectorF* a, const AVXVectorF* b) {
  if (result->_size != a->_size || result->_size != b->_size) {
    throw std::runtime_error("Not the same size in max");
  }
  std::size_t iterations = a->getIterations();
  for (std::size_t i = 0; i < iterations; i++) {
    AVL_FloatVectorType a_data = AVL_load_ps (a->_data + i * AVL_VectorSize);
    AVL_FloatVectorType b_data = AVL_load_ps (b->_data + i * AVL_VectorSize);
    AVL_MaskType result_data = AVL_cmp_ps_mask(a_data, b_data, _CMP_LE_OQ);
    AVL_store_mask(result->_data + i, result_data);

  }
}





}
