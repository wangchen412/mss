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

# set(mss mss_core mss_input)

function(mss_add_test test_name)
  add_executable(test_${test_name} test_${test_name}.cc)
  target_link_libraries(test_${test_name} gtest gtest_main complex_bessel)
  add_test(NAME test_${test_name} COMMAND test_${test_name})
endfunction()
