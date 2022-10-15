#Currently the lastest CMake version
cmake_minimum_required (VERSION 3.9)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)




macro (create_test target)
	ADD_EXECUTABLE( ${target} ${target}.cxx )
	target_link_libraries(${target}   AVX512VectorLib  ${CMAKE_SHARED_LINKER_FLAGS}  )
	set_target_properties(${target} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/tests)
  	add_test (NAME ${target} COMMAND $<TARGET_FILE:${target}>)
endmacro (create_test)

create_test(atan)
