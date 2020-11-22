add_subdirectory(io)
add_subdirectory(geometry/octree)
add_executable(nebula_cpu_dynamic cpu_dynamic.cpp
		$<TARGET_OBJECTS:io>)
#target_link_libraries(nebula_cpu_mt ${LIBRARIES} Threads::Threads)
#target_link_libraries(nebula_cpu_edep ${LIBRARIES} Threads::Threads)
target_link_libraries(nebula_cpu_dynamic ${LIBRARIES} Threads::Threads)
