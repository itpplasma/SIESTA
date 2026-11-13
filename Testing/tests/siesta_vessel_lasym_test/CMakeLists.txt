#  Copy input files from the source to the build directory.
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/siesta.jcf
                ${CMAKE_CURRENT_BINARY_DIR}/siesta.jcf
                @ONLY)
foreach (file input.test.vmec vessel.test)
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/${file}
                    ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY)
endforeach ()
foreach (file coils.test)
    configure_file (${CMAKE_CURRENT_SOURCE_DIR}/../siesta_vessel_test/${file}
                    ${CMAKE_CURRENT_BINARY_DIR}/${file}
                    COPYONLY)
endforeach ()

#  Fixed-boundary VMEC preserves the synthetic asymmetric boundary modes.
add_test (NAME    siesta_vessel_lasym_mgrid_test
          COMMAND $<TARGET_PROPERTY:mgrid,BINARY_DIR>/mgrid coils.test)
set_tests_properties (siesta_vessel_lasym_mgrid_test
                      PROPERTIES FIXTURES_SETUP siesta_vessel_lasym_mgrid)
add_test (NAME    siesta_vessel_lasym_vmec_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xvmec,BINARY_DIR>/xvmec input.test.vmec)
set_tests_properties (siesta_vessel_lasym_vmec_test
                      PROPERTIES
                      DEPENDS siesta_vessel_lasym_mgrid_test
                      FIXTURES_REQUIRED siesta_vessel_lasym_mgrid
                      FIXTURES_SETUP siesta_vessel_lasym_vmec
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>)
add_test (NAME    siesta_vessel_lasym_siesta_test
          COMMAND $<TARGET_PROPERTY:stell,MPIEXEC_EXECUTABLE> $<TARGET_PROPERTY:stell,MPIEXEC_NUMPROC_FLAG> $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS> $<TARGET_PROPERTY:xsiesta,BINARY_DIR>/xsiesta)
set_tests_properties (siesta_vessel_lasym_siesta_test
                      PROPERTIES
                      DEPENDS siesta_vessel_lasym_vmec_test
                      FIXTURES_REQUIRED siesta_vessel_lasym_vmec
                      PROCESSORS $<TARGET_PROPERTY:stell,MPIEXEC_MAX_NUMPROCS>)
