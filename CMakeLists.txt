cmake_minimum_required(VERSION 3.27)
project(ASTRA C)

include_directories(/usr/local/opt/curl/include/)

set(CMAKE_C_STANDARD 11)

set(CURL_LIBRARY "-lcurl")
find_package(CURL REQUIRED)

add_executable(ASTRA main.c)

include_directories(${CURL_INCLUDE_DIR})
target_link_libraries(ASTRA2 ${CURL_LIBRARIES})