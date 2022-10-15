#ifndef _AVX512VECTORLIB_H_INCLUDED
#define _AVX512VECTORLIB_H_INCLUDED

//system includes
#include <immintrin.h>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <cstring>
#include <iomanip>
#include <limits>
//library includes

//own includes


namespace InstabilityDetectionLib {


#ifdef __AVX512F__
	#define AVL_USE_AVX512
	//#warning "Using AVX512"
	//#pragma message("Using AVX512")
#endif

#ifndef __AVX512F__
#ifdef __AVX2__
	#define AVL_USE_AVX2
	//#pragma message("Using AVX2")
	//#warning "Using AVX2"
#else
	#warning "No AVX detected"
#endif
#endif

#ifdef AVL_USE_AVX512

#define AVL_FloatVectorType __m512
#define AVL_IntegerVectorType __m512i
#define AVL_IntegerVectorConversionType __m256i
#define AVL_MaskType __mmask16
#define AVL_VectorSize 16
#define AVL_Alignment 64

#define AVL_load_ps _mm512_load_ps
#define AVL_cmp_ps_mask _mm512_cmp_ps_mask
#define AVL_store_mask _store_mask16
#define AVL_set1_ps _mm512_set1_ps
#define AVL_max_ps _mm512_max_ps
#define AVL_min_ps _mm512_min_ps
#define AVL_store_ps _mm512_store_ps
#define AVL_andnot_ps _mm512_andnot_ps
#define AVL_kor _mm512_kor
#define AVL_mask_blend_ps _mm512_mask_blend_ps
#define AVL_setzero_ps _mm512_setzero_ps
#define AVL_mask_max_ps _mm512_mask_max_ps
#define AVL_mask_min_ps _mm512_mask_min_ps
#define AVL_kand _mm512_kand
#define AVL_mask_add_ps _mm512_mask_add_ps
#define AVL_add_ps _mm512_add_ps
#define AVL_mask_mov_ps _mm512_mask_mov_ps
#define AVL_mask_sub_ps _mm512_mask_sub_ps
#define AVL_fmadd_ps _mm512_fmadd_ps
#define AVL_mul_ps _mm512_mul_ps
#define AVL_fmsub_ps _mm512_fmsub_ps
#define AVL_knot _mm512_knot
#define AVL_mask_div_ps _mm512_mask_div_ps
#define AVL_sub_ps _mm512_sub_ps
#define AVL_rcp14_ps _mm512_rcp14_ps
#define AVL_div_ps _mm512_div_ps
#define AVL_mask_store_ps _mm512_mask_store_ps
#define AVL_cvtepi16_epi32 _mm512_cvtepi16_epi32
#define AVL_maskz_loadu_epi16 _mm256_maskz_loadu_epi16
#define AVL_load_mask _load_mask16
#define AVL_castsi_ps _mm512_castsi512_ps
#define AVL_xor_si _mm512_xor_si512
#define AVL_castps_si _mm512_castps_si512
#define AVL_and_si _mm512_and_si512
#define AVL_andnot_si _mm512_andnot_si512
#define AVL_srli_epi32 _mm512_srli_epi32
#define AVL_sub_epi32 _mm512_sub_epi32
#define AVL_mask_fmsub_ps _mm512_mask_fmsub_ps
#define AVL_fnmadd_ps _mm512_fnmadd_ps
#define AVL_or_si _mm512_or_si512
#define AVL_cvtepi32_ps _mm512_cvtepi32_ps
#define AVL_floor_ps _mm512_floor_ps
#define AVL_cvttps_epi32 _mm512_cvttps_epi32
#define AVL_add_epi32 _mm512_add_epi32
#define AVL_slli_epi32 _mm512_slli_epi32
#define AVL_cmpeq_epi32_mask _mm512_cmpeq_epi32_mask
#define AVL_mask_blend_epi32 _mm512_mask_blend_epi32
#define AVL_cvtu32_mask _cvtu32_mask16
#define AVL_maskz_rsqrt14_ps _mm512_maskz_rsqrt14_ps
#define AVL_abs_ps	_mm512_abs_ps


/* declare some AVX512 constants -- why can't I figure a better way to do that? */
#define _PS512_CONST(Name, Val)                                            \
  static const ALIGN32_BEG float _ps512_##Name[AVL_VectorSize] ALIGN32_END = { Val, Val, Val, Val, Val, Val, Val, Val, Val, Val, Val, Val, Val, Val, Val, Val }
#define _PI32_CONST512(Name, Val)                                            \
  static const ALIGN32_BEG int _pi32_512_##Name[AVL_VectorSize] ALIGN32_END = { Val, Val, Val, Val, Val, Val, Val, Val, Val, Val, Val, Val, Val, Val, Val, Val }
#define _PS512_CONST_TYPE(Name, Type, Val)                                 \
  static const ALIGN32_BEG Type _ps512_##Name[AVL_VectorSize] ALIGN32_END = { Val, Val, Val, Val, Val, Val, Val, Val, Val, Val, Val, Val, Val, Val, Val, Val }

#endif

#ifdef AVL_USE_AVX2
#define AVL_FloatVectorType __m256
#define AVL_IntegerVectorType __m256i
#define AVL_IntegerVectorConversionType __m128i
#define AVL_MaskType __mmask8
#define AVL_VectorSize 8
#define AVL_Alignment 32

#define AVL_load_ps _mm256_load_ps
#define AVL_cmp_ps_mask _mm256_cmp_ps_mask
#define AVL_store_mask _store_mask8
#define AVL_set1_ps _mm256_set1_ps
#define AVL_max_ps _mm256_max_ps
#define AVL_min_ps _mm256_min_ps
#define AVL_store_ps _mm256_store_ps
#define AVL_andnot_ps _mm256_andnot_ps
#define AVL_kor _kor_mask8
#define AVL_mask_blend_ps _mm256_mask_blend_ps
#define AVL_setzero_ps _mm256_setzero_ps
#define AVL_mask_max_ps _mm256_mask_max_ps
#define AVL_mask_min_ps _mm256_mask_min_ps
#define AVL_kand _kand_mask8
#define AVL_mask_add_ps _mm256_mask_add_ps
#define AVL_add_ps _mm256_add_ps
#define AVL_mask_mov_ps _mm256_mask_mov_ps
#define AVL_mask_sub_ps _mm256_mask_sub_ps
#define AVL_fmadd_ps _mm256_fmadd_ps
#define AVL_mul_ps _mm256_mul_ps
#define AVL_fmsub_ps _mm256_fmsub_ps
#define AVL_knot _knot_mask8
#define AVL_mask_div_ps _mm256_mask_div_ps
#define AVL_sub_ps _mm256_sub_ps
#define AVL_rcp14_ps _mm256_rcp14_ps
#define AVL_div_ps _mm256_div_ps
#define AVL_mask_store_ps _mm256_mask_store_ps
#define AVL_cvtepi16_epi32 _mm256_cvtepi16_epi32
#define AVL_maskz_loadu_epi16 _mm_maskz_loadu_epi16
#define AVL_load_mask _load_mask8
#define AVL_castsi_ps _mm256_castsi256_ps
#define AVL_xor_si _mm256_xor_si256
#define AVL_castps_si _mm256_castps_si256
#define AVL_and_si _mm256_and_si256
#define AVL_andnot_si _mm256_andnot_si256
#define AVL_srli_epi32 _mm256_srli_epi32
#define AVL_sub_epi32 _mm256_sub_epi32
#define AVL_mask_fmsub_ps _mm256_mask_fmsub_ps
#define AVL_fnmadd_ps _mm256_fnmadd_ps
#define AVL_or_si _mm256_or_si256
#define AVL_cvtepi32_ps _mm256_cvtepi32_ps
#define AVL_floor_ps _mm256_floor_ps
#define AVL_cvttps_epi32 _mm256_cvttps_epi32
#define AVL_add_epi32 _mm256_add_epi32
#define AVL_slli_epi32 _mm256_slli_epi32
#define AVL_cmpeq_epi32_mask _mm256_cmpeq_epi32_mask
#define AVL_mask_blend_epi32 _mm256_mask_blend_epi32
#define AVL_cvtu32_mask _cvtu32_mask8
#define AVL_maskz_rsqrt14_ps _mm256_maskz_rsqrt14_ps
#define AVL_abs_ps	_mm256_abs_ps

inline __m256 _mm256_abs_ps(__m256 x) {
    static const __m256 sign_mask = _mm256_set1_ps(-0.f); // -0.f = 1 << 31
    return _mm256_andnot_ps(sign_mask, x);
}


/* declare some AVX512 constants -- why can't I figure a better way to do that? */
#define _PS512_CONST(Name, Val)                                            \
  static const ALIGN32_BEG float _ps512_##Name[AVL_VectorSize] ALIGN32_END = { Val, Val, Val, Val, Val, Val, Val, Val  }
#define _PI32_CONST512(Name, Val)                                            \
  static const ALIGN32_BEG int _pi32_512_##Name[AVL_VectorSize] ALIGN32_END = { Val, Val, Val, Val, Val, Val, Val, Val }
#define _PS512_CONST_TYPE(Name, Type, Val)                                 \
  static const ALIGN32_BEG Type _ps512_##Name[AVL_VectorSize] ALIGN32_END = { Val, Val, Val, Val, Val, Val, Val, Val }


#endif





//typedef AVL_FloatVectorType AVL_FloatVectorType; // vector of 8 float (avx)

// prototypes
inline AVL_FloatVectorType log512_ps(AVL_FloatVectorType x);
inline AVL_FloatVectorType exp512_ps(AVL_FloatVectorType x);
inline AVL_FloatVectorType sin512_ps(AVL_FloatVectorType x);
inline AVL_FloatVectorType cos512_ps(AVL_FloatVectorType x);
inline void sincos512_ps(AVL_FloatVectorType x, AVL_FloatVectorType *s, AVL_FloatVectorType *c);

int roundUp(std::size_t numToRound, std::size_t multiple);



class AVXVectorMask {
public:
	AVXVectorMask(const std::size_t& size);
	~AVXVectorMask();
	bool operator [](int i) const {
		AVL_MaskType result_data = _data[i / 16];
		return result_data & 1 << (i % 16);

	}
	void print() {
		for (std::size_t i = 0; i < 8; i++) {
			std::cout << (*this)[i] << ",";
		}
		std::cout << std::endl;
	}
	void setData(const AVXVectorMask* data);
	void setData(const bool& data);
	//result=a|b
	static  void OR(const AVXVectorMask* result, const AVXVectorMask* a, const AVXVectorMask* b);
	static  void AND(const AVXVectorMask* result, const AVXVectorMask* a, const AVXVectorMask* b);
	//result = not a
	static  void NOT(const AVXVectorMask* result, const AVXVectorMask* a);

protected:
private:
	std::size_t getIterations() const;
	std::size_t _size = 0;
	AVL_MaskType* _data;
	friend class AVXVectorF;
};


class AVXVectorF {
public:
	AVXVectorF(const std::size_t& size);
	~AVXVectorF();
	float operator [](int i) const    {return _data[i];}
	float & operator [](int i) {return _data[i];}

	const float* getPtr(){
		return _data;
	}

	void setData(const float& data);

	void setData(const float* data);

	void setData(const int16_t* data);

	void setData(const std::vector<float>& data);

	void setData(const AVXVectorF* data);
	void setData(const AVXVectorMask* mask, const float& a);

	void print();

	//result=a+b
	static  void add(const AVXVectorF* result, const AVXVectorF* a, const AVXVectorF* b);
	//result=a+b
	static  void add(const AVXVectorF* result, const AVXVectorF* a, const float& b);
	//if mask:
	//	result=a+b
	//else
	//result=result
	static  void add_mask(const AVXVectorF* result, const AVXVectorMask* mask, const AVXVectorF* a, const float& b);
	//result=a-b
	static  void sub(const AVXVectorF* result, const AVXVectorF* a, const AVXVectorF* b);
	//result=a*b
	static  void mult(const AVXVectorF* result, const AVXVectorF* a, const AVXVectorF* b);
	//result=a*b
	static  void mult(const AVXVectorF* result, const AVXVectorF* a, const float& b);
	//result=a/b
	static  void div(const AVXVectorF* result, const AVXVectorF* a, const AVXVectorF* b);
	//result=a/b
	static  void div(const AVXVectorF* result, const AVXVectorF* a, const float& b);
	//result=a/b  x/0=inf -x/0=-inf
	static  void divWithInf(const AVXVectorF* result, const AVXVectorF* a, const AVXVectorF* b);
	// result=a*b+c
	static  void fmadd(const AVXVectorF* result, const AVXVectorF* a, const AVXVectorF* b, const AVXVectorF* c);
	// result=a*b+c
	static  void fmadd(const AVXVectorF* result, const AVXVectorF* a, const float& b, const AVXVectorF* c);
// result=a*b-c
	static  void fmsub(const AVXVectorF* result, const AVXVectorF* a, const AVXVectorF* b, const AVXVectorF* c);
	//result=ln(a)
	static  void ln(const AVXVectorF* result, const AVXVectorF* a);
	//result=exp(a)
	static  void exp(const AVXVectorF* result, const AVXVectorF* a);
	//result=sqrt(a)
	static  void sqrt(const AVXVectorF* result, const AVXVectorF* a);

	static void complexDiv(const AVXVectorF* resultReal, const AVXVectorF* resultImag, const AVXVectorF* a, const AVXVectorF* b, const AVXVectorF* c, const AVXVectorF* d);
	//result=atan2(a,b)
	static void atan2(const AVXVectorF * result, const AVXVectorF * a, const AVXVectorF * b);
	//returns the number of occurences of b in a
	static  std::size_t count(const AVXVectorF * a,const float& b);
	//result=sum(a)
	static  float sum(const AVXVectorF * a);
	//result=min(a)
	static  float min(const AVXVectorF * a);
	//result=min(a) where a is not zero
	static  float minNotZero(const AVXVectorF * a);
	//result=max(a)
	static  float max(const AVXVectorF * a);
	//removes any inf and nan
	static  void rmInfNan(const AVXVectorF * a);
	//result=min(a,b)
	static  void min(const AVXVectorF * result, const AVXVectorF * a, const AVXVectorF * b);
	//result=max(a,b)
	static  void max(const AVXVectorF * result, const AVXVectorF * a, const AVXVectorF * b);
	// result=a==b
	static  void equal(const AVXVectorMask * result, const AVXVectorF * a, const AVXVectorF * b);
	// result=a==b
	static  void equal(const AVXVectorMask * result, const AVXVectorF * a, const float & b);
	//result= a>b
	static  void greater_than(const AVXVectorMask * result, const AVXVectorF * a, const AVXVectorF * b);
	//result= a>b
	static  void greater_than(const AVXVectorMask * result, const AVXVectorF * a, const float & b);
	//result= a>=b
	static  void greater_equal_than(const AVXVectorMask * result, const AVXVectorF * a, const AVXVectorF * b);
	//result= a<b
	static  void less_than(const AVXVectorMask * result, const AVXVectorF * a, const AVXVectorF * b);
	//result= a<=b
	static  void less_equal_than(const AVXVectorMask * result, const AVXVectorF * a, const AVXVectorF * b);

	static void print(AVL_FloatVectorType & data) {
		float result_data[AVL_VectorSize];
		AVL_store_ps(result_data , data);
		for (std::size_t i = 0; i < AVL_VectorSize; i++) {
			std::cout << result_data[i] << " ";
		}
		std::cout << std::endl;
	}

	static void print(AVL_MaskType & data) {
		for (std::size_t i = 0; i < AVL_VectorSize; i++) {
			std::cout << ((data & 1 << i) != 0) << ", ";
		}
		std::cout << std::endl;
	}


protected:
private:
	AVL_MaskType getLastMask() const;
	//required for sum,min,max where we should not use all the values in the last vector
	std::size_t getIterations() const;

	std::size_t _size = 0;
	float* _data;
};


}

#endif
