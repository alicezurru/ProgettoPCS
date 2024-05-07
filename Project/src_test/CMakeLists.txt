list (APPEND src_sources "${CMAKE_CURRENT_SOURCE_DIR}/Fractures.cpp")
list (APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/Fractures.hpp")

list (APPEND src_sources "${CMAKE_CURRENT_SOURCE_DIR}/Utils.cpp")
list (APPEND src_headers "${CMAKE_CURRENT_SOURCE_DIR}/Utils.hpp")

list (APPEND src_test_headers "${CMAKE_CURRENT_SOURCE_DIR}/Test.hpp")

list (APPEND src_test_includes "${CMAKE_CURRENT_SOURCE_DIR}")

set (src_test_sources ${src_test_sources} PARENT_SCOPE)
set (src_test_headers ${src_test_headers} PARENT_SCOPE)
set (src_test_includes ${src_test_includes} PARENT_SCOPE)
