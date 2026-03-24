#!/bin/bash

# Copyright 2023 University of Warwick
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# Laura++ package authors:
# John Back
# Paul Harrison
# Thomas Latham

# apply clang-format to all source files
find . \( -name "*.hh" -o -name "*.cc" -o -name "*.C" \) -exec clang-format -i '{}' \;

# fix comment sections in all source files (clang-format does not touch these)
find . \( -name "*.hh" -o -name "*.cc" -o -name "*.C" \) -exec bash fix-comments.sh '{}' \;

# remove tabs from CMake files
find . \( -name CMakeLists.txt -o -name "*.cmake" \) -exec sed -i -e 's/\t/    /g' '{}' \;

# remove trailing whitespace in all Markdown files
find . -name "*.md" -exec sed -i -e 's/[[:space:]]\+$//' '{}' \;

