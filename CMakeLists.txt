# ----------------------------------------------------------------------
#
# Copyright © 2017 mss authors.
#
# mss is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License,
# or (at your option) any later version.
#
# mss is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# ----------------------------------------------------------------------

cmake_minimum_required(VERSION 3.7)

set(CMAKE_CXX_COMPILER icpc)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# find_package(Eigen3 REQUIRED)
# include_directories(${EIGEN3_INCLUDE_DIR})
add_definitions(-D EIGEN_USE_MKL_ALL)

include_directories(/home/wangchen/.third_party)
link_directories(/home/wangchen/.third_party)

project(mss)

# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

set(DFLAGS "-Wall -mkl=sequential -fp-model precise")
set(RFLAGS "-Wall -mkl=parallel -fopenmp -parallel -fp-model precise")

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${DFLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${RFLAGS}")

enable_testing()

add_subdirectory(com)
add_subdirectory(tests)
