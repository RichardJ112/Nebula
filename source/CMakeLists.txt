add_subdirectory(io)
add_subdirectory(geometry/octree)


#add_executable(nebula_cpu_mt cpu_mt_main.cpp
#		$<TARGET_OBJECTS:io>
#		$<TARGET_OBJECTS:octree>)
#add_executable(nebula_cpu_edep cpu_energydep.cpp
#		$<TARGET_OBJECTS:io>
#		$<TARGET_OBJECTS:octree>)
add_executable(nebula_cpu_custom cpu_custom.cpp
		$<TARGET_OBJECTS:io>)
#target_link_libraries(nebula_cpu_mt ${LIBRARIES} Threads::Threads)
#target_link_libraries(nebula_cpu_edep ${LIBRARIES} Threads::Threads)
target_link_libraries(nebula_cpu_custom ${LIBRARIES} Threads::Threads)
