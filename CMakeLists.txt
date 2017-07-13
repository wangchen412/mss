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
#set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 14)


project(mss)

enable_testing()

include_directories("~/Programming/eigen/")
include_directories("~/.third_party/")
link_directories("~/.third_party/")
add_subdirectory(src)
add_subdirectory(tests)
