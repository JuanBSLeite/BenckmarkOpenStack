#!/usr/bin/env python

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

import sys
import re

if len(sys.argv) != 2 :
    print(f'Usage: {sys.argv[0]} <filename>')
    sys.exit(1)

filename = sys.argv[1]

lines = []
with open(filename) as checkfile :
    lines = checkfile.readlines()

current_opening = None
comment_blocks = []

# find all multi-line blocks of Doxygen comments, delimited by /*! and */
for i in range( 0, len(lines) ) :

    myline = lines[i]

    opening_match = re.search( '^\s*[/][*][!]', myline )

    if opening_match :
        if current_opening :
            raise Exception( f'Found a second opening to a comment block on line {i+1} while still looking for the closing of the block opened on line {current_opening+1}' )
        current_opening = i

    if not current_opening :
        continue

    closing_match = re.search( '^\s*[*][/]', myline )

    if closing_match :
        comment_blocks.append( (current_opening, i) )
        current_opening = None

if len(comment_blocks) == 0 :
    print('Found no comment blocks!')
    sys.exit(0)

need_new_line = []
changesMade = False
for first, final in comment_blocks :

    first_line = lines[first]
    final_line = lines[final]

    first_indentation_level = first_line.index('/')
    final_indentation_level = final_line.index('*')

    # if necessary, fix the opening line
    if first_line.endswith(' ') :
        lines[first] = first_line.rstrip() + '\n'
        changesMade = True

    # if necessary, fix the closing line
    if final_indentation_level != first_indentation_level or final_line.endswith(' ') :
        newline = first_indentation_level*' ' + final_line.strip()
        lines[final] = newline + '\n'
        changesMade = True

    block_indentation_level = first_indentation_level + 4

    # check and, if necessary, fix all intervening lines
    for i in range( first+1, final ) :

        myline = lines[i].rstrip('\n')
        if len(myline) == 0 :
            continue

        stripped_line = myline.lstrip()
        if len(stripped_line) == 0 :
            lines[i] = '\n'
            changesMade = True
            continue

        indentation_level = len(myline) - len(stripped_line)

        if indentation_level != block_indentation_level or stripped_line.endswith(' ') :
            newline = block_indentation_level*' ' + stripped_line.rstrip()
            lines[i] = newline + '\n'
            changesMade = True

        # identify where we need to insert blank lines between description and parameters (and similar)
        if i+1 < final :
            nextline = lines[i+1]
            if not stripped_line.startswith('\\') and nextline.lstrip().startswith('\\') :
                need_new_line.append(i+1)
                changesMade = True

# if no changes mode (or required in the form of extra blank lines) we can exit
if not changesMade :
    sys.exit(0)

# insert any required blank lines
need_new_line.sort(reverse=True)
for i in need_new_line :
    lines.insert(i, '\n')

# write out modified file
with open(filename+'.new', 'w') as outputfile :
    outputfile.writelines(lines)
