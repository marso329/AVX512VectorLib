#include "AVX512VectorLib.hpp"
#include <cmath>


using namespace InstabilityDetectionLib;


int main() {
	std::size_t length=1024;

	AVXVectorF* x = new AVXVectorF(length);
	AVXVectorF* y = new AVXVectorF(length);
	AVXVectorF* result = new AVXVectorF(length);
	for(std::size_t i=0;i<length;i++){
		(*x)[i]=4.0- 2.0f*static_cast<float>(i)/static_cast<float>(length);
		(*y)[i]=-1.0+ 2.0f*static_cast<float>(i)/static_cast<float>(length);
	}
	AVXVectorF::atan2(result,y,x);
	float max_error=0;
	for(std::size_t i=0;i<length;i++){
		float own_result=(*result)[i];
		float std_result=atan2((*y)[i],(*x)[i]);
		float error=abs(std_result-own_result );
		std::cout<<"Own result: "<<own_result<<" std_result: "<<std_result<<" error: "<<error<<" index: "<<i<<" y: "<<(*y)[i]<<" x: "<<(*x)[i]<<std::endl;
		max_error=std::max(max_error,error);

	}
	std::cout<<max_error<<std::endl;

}
