cmake_minimum_required(VERSION 3.10)
project(AneurysmPlugin)

set(CMAKE_CXX_STANDARD 14)

# Force the architecture to match your FEBioStudio (x86_64)
set(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE STRING "" FORCE)

add_library(AneurysmPlugin SHARED AneurysmPlugin.cpp)

# Use the SDK path you verified
set(FEBIO_DIR "/Applications/FEBioStudio/sdk" CACHE PATH "Path to FEBio SDK" FORCE)

target_include_directories(AneurysmPlugin PRIVATE "${FEBIO_DIR}/include")
target_link_directories(AneurysmPlugin PRIVATE "${FEBIO_DIR}/lib") 
target_link_libraries(AneurysmPlugin PRIVATE fecore febiomech)

# --- THE MAC FIX ---
# This tells the plugin to look for FEBio libraries in the same folder as the plugin
# or in the standard FEBio library locations.
set_target_properties(AneurysmPlugin PROPERTIES 
    PREFIX "" 
    SUFFIX ".dylib"
    INSTALL_RPATH "@loader_path;@loader_path/../lib;${FEBIO_DIR}/lib"
    BUILD_WITH_INSTALL_RPATH TRUE
)