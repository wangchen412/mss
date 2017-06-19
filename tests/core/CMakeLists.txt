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

add_executable(test_math test_math.cc)
target_link_libraries(test_math mss_core gtest gtest_main)
add_test(NAME test_math COMMAND test_math)

add_executable(test_cs test_cs.cc)
target_link_libraries(test_cs mss_core gtest gtest_main)
add_test(NAME test_cs COMMAND test_cs)

add_executable(test_state test_state.cc)
target_link_libraries(test_state mss_core gtest gtest_main)
add_test(NAME test_state COMMAND test_cs)
