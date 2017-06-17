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

project(mss)
set(CMAKE_CXX_STANDARD 14)

enable_testing()

include_directories("~/Programming/eigen/")
add_subdirectory(src)
add_subdirectory(tests)


# add_custom_target(test)
# macro(run_test test_target)
#   add_custom_target(${test_target}_runtest
#       COMMAND ${test_target}
#       DEPENDS ${test_target}
#       WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
#   add_dependencies(test ${test_target}_runtest)
# endmacro()
