find_package(Doxygen)

if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plll.sty.in ${CMAKE_CURRENT_BINARY_DIR}/html/plll.sty @ONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/plll.sty.in ${CMAKE_CURRENT_BINARY_DIR}/latex/plll.sty @ONLY)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/images/makefile.in ${CMAKE_CURRENT_BINARY_DIR}/images/makefile @ONLY)
    add_custom_target(docimages # ALL # uncomment if documentation should be made automatically
                      COMMAND make
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/images
                      COMMENT "Generating images for documentation" VERBATIM)
    add_custom_target(doc # ALL # uncomment if documentation should be made automatically
                      COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                      DEPENDS docimages
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      COMMENT "Generating documentation (with Doxygen)" VERBATIM)
else(DOXYGEN_FOUND)
    message(STATUS "Doxygen was not found, so documentation will not be created.")
endif(DOXYGEN_FOUND)
