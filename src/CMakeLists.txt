macro(add_simple_observer observer_name)
  add_observer(${observer_name} "${observer_name}.cpp" "${CMAKE_SOURCE_DIR}/include/${PROJECT_NAME}/${observer_name}.h")
  # mc_observers doesn't actually depend on mc_control, but the observer's implementation does
  target_link_libraries(${observer_name} PUBLIC mc_rtc::mc_control)
endmacro()

macro(add_so_observer observer_name)
  add_simple_observer(${observer_name})
  target_link_libraries(${observer_name} PUBLIC state-observation::state-observation)
endmacro()

macro(add_observer_with_filter observer_name)
  add_observer(${observer_name} "${observer_name}.cpp;filtering.cpp" "${CMAKE_SOURCE_DIR}/include/${PROJECT_NAME}/${observer_name}.h;${CMAKE_SOURCE_DIR}/include/${PROJECT_NAME}/filtering.h")
  target_link_libraries(${observer_name} PUBLIC mc_rtc::mc_control)
endmacro()

add_so_observer(ContactEstimation)


install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/${PROJECT_NAME} DESTINATION include)
