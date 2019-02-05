include_directories(${PROJECT_SOURCE_DIR}/plll/include ${PROJECT_BINARY_DIR}/plll/include ${PROJECT_SOURCE_DIR}/plll/include-internal)

set(LIBPLLL_SOURCES src/arithmetic/arithmetic.cpp
                    src/arithmetic/arithmetic-gmp.cpp
                    src/arithmetic/rational.cpp
                    src/arithmetic/primes.cpp
                    src/arithmetic/factor.cpp
                    src/arithmetic/ddqd-wrapper.cpp
                    src/linalg/linalg.cpp
                    src/lattices/lll2.cpp
                    src/lattices/lll2-multiplexer.cpp
                    src/lattices/transform-impl.cpp
                    src/lattices/compilewrappers/lll2-compunit-bigint-rational-gs.cpp
                    src/lattices/compilewrappers/lll2-compunit-bigint-real-gs.cpp
                    src/lattices/compilewrappers/lll2-compunit-bigint-longdouble-gs.cpp
                    src/lattices/compilewrappers/lll2-compunit-bigint-double-gs.cpp
                    src/lattices/compilewrappers/lll2-compunit-bigint-doubledouble-gs.cpp
                    src/lattices/compilewrappers/lll2-compunit-bigint-quaddouble-gs.cpp
                    src/lattices/compilewrappers/lll2-compunit-bigint-rational-red.cpp
                    src/lattices/compilewrappers/lll2-compunit-bigint-real-red.cpp
                    src/lattices/compilewrappers/lll2-compunit-bigint-longdouble-red.cpp
                    src/lattices/compilewrappers/lll2-compunit-bigint-double-red.cpp
                    src/lattices/compilewrappers/lll2-compunit-bigint-doubledouble-red.cpp
                    src/lattices/compilewrappers/lll2-compunit-bigint-quaddouble-red.cpp
                    src/lattices/compilewrappers/lll2-compunit-bigint-rational-enum.cpp
                    src/lattices/compilewrappers/lll2-compunit-bigint-real-enum.cpp
                    src/lattices/compilewrappers/lll2-compunit-bigint-longdouble-enum.cpp
                    src/lattices/compilewrappers/lll2-compunit-bigint-double-enum.cpp
                    src/lattices/compilewrappers/lll2-compunit-bigint-doubledouble-enum.cpp
                    src/lattices/compilewrappers/lll2-compunit-bigint-quaddouble-enum.cpp
                    src/lattices/compilewrappers/lll2-compunit-longint-rational-gs.cpp
                    src/lattices/compilewrappers/lll2-compunit-longint-real-gs.cpp
                    src/lattices/compilewrappers/lll2-compunit-longint-longdouble-gs.cpp
                    src/lattices/compilewrappers/lll2-compunit-longint-double-gs.cpp
                    src/lattices/compilewrappers/lll2-compunit-longint-doubledouble-gs.cpp
                    src/lattices/compilewrappers/lll2-compunit-longint-quaddouble-gs.cpp
                    src/lattices/compilewrappers/lll2-compunit-longint-rational-red.cpp
                    src/lattices/compilewrappers/lll2-compunit-longint-real-red.cpp
                    src/lattices/compilewrappers/lll2-compunit-longint-longdouble-red.cpp
                    src/lattices/compilewrappers/lll2-compunit-longint-double-red.cpp
                    src/lattices/compilewrappers/lll2-compunit-longint-doubledouble-red.cpp
                    src/lattices/compilewrappers/lll2-compunit-longint-quaddouble-red.cpp
                    src/lattices/compilewrappers/lll2-compunit-longint-rational-enum.cpp
                    src/lattices/compilewrappers/lll2-compunit-longint-real-enum.cpp
                    src/lattices/compilewrappers/lll2-compunit-longint-longdouble-enum.cpp
                    src/lattices/compilewrappers/lll2-compunit-longint-double-enum.cpp
                    src/lattices/compilewrappers/lll2-compunit-longint-doubledouble-enum.cpp
                    src/lattices/compilewrappers/lll2-compunit-longint-quaddouble-enum.cpp
                    src/utilities/arguments.cpp
                    src/utilities/timer.cpp
                    src/utilities/keccak.cpp
                    src/utilities/dlalloc.cpp
                    src/utilities/myalloc.cpp
                    src/utilities/taskmanager.cpp)

add_library(libplll-static STATIC ${LIBPLLL_SOURCES})
set_target_properties(libplll-static PROPERTIES OUTPUT_NAME plll)
set_target_properties(libplll-static PROPERTIES CLEAN_DIRECT_OUTPUT 1)
target_link_libraries(libplll-static ${libraries})

add_library(libplll-shared SHARED ${LIBPLLL_SOURCES})
set_target_properties(libplll-shared PROPERTIES OUTPUT_NAME plll)
set_target_properties(libplll-shared PROPERTIES CLEAN_DIRECT_OUTPUT 1)
target_link_libraries(libplll-shared ${libraries_shared})

install(TARGETS libplll-static libplll-shared DESTINATION lib)
install(FILES include/plll.hpp DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/plll/include/plll/config.hpp
              include/plll/arguments.hpp
              include/plll/arithmetic.hpp
              include/plll/arithmetic-expressions.hpp
              include/plll/arithmetic-gmp.hpp
              include/plll/arithmetic-gmp-conv.hpp
              include/plll/arithmetic-gmp-iops.hpp
              include/plll/arithmetic-gmp-rops.hpp
              include/plll/arithmetic-nint.hpp
              include/plll/arithmetic-nint-conv.hpp
              include/plll/documentation.hpp
              include/plll/helper.hpp
              include/plll/linalg.hpp
              include/plll/matrix.hpp
              include/plll/matrix-mem.hpp
              include/plll/matrix-ops.hpp
              include/plll/matrix-ops2.hpp
              include/plll/rational.hpp
              include/plll/rational-conv.hpp
              include/plll/rational-ops.hpp
              DESTINATION include/plll)

