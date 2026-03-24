
list(APPEND TEST_SOURCES
    TestCovariant
    TestCovariant2
    TestNewKinematicsMethods
    )

foreach( _test ${TEST_SOURCES})
    add_executable(${_test} ${_test}.cc)
    target_link_libraries(${_test} PRIVATE Laura++)
    #install(TARGETS ${_test} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endforeach()

