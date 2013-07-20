# author   : Johann-Mattis List
# email    : mattis.list@uni-marburg.de
# created  : 2013-07-17 10:40
# modified : 2013-07-17 10:40
"""
Module handels all global parameters used in a LingPy session.
"""

__author__="Johann-Mattis List"
__date__="2013-07-17"

# builtin imports
from datetime import datetime,date
import os

# set the absolute path to lingpy
_abs_path = os.path.split(
        os.path.abspath(
            __file__
            )
        )[0]

# initialize rcParams with filename and timestamp
rcParams = dict(
        _path = _abs_path,
        filename                   = 'lingpy-'+str(date.today()),
        timestamp                  = str(datetime.today()),
        answer_yes = ['y','Y','j','J','yes'],
        verbose = False,
        debug = False
        )

# messages
messages = dict(
        M_file_written = "Data has been written to file <{0}>.",
        M_no_errors_in_data = "No obvious errors found in the data.",
        M_random_alignments = "Calculating random alignments for pair {0} / {1}.",
        M_alignments = "Calculating alignments for pair {0} / {1}.",
        )
# adjust the prefix for general-purpose messages
for k,v in messages.items():
    messages[k] = "[i] "+v
# update rcParams
rcParams.update(messages)

# these are general warnings, all prefixed by "warning" in rcParams
warnings = dict(
        W_empty_cons       = "There are emtpy segments in the consensus.",
        W_failed_cons      = "Failed to compute the consensus string.",
        W_deprecation      = "Use of '{0}' is deprecated, use '{1}' instead.",
        W_missing_module   = "Module '{0}' could not be loaded. Some methods may not work properly.",
        W_identical_scorer = "An identical scoring function has already been calculated, force recalculation by setting 'force' to 'True'.",
        W_overwrite_scorer = "A different scoring function has already been calculated, overwriting previous settings.",
        W_zero_division    = "Zero-division error encountered in '{0}' and '{1}'."
        )
# adjust the prefix for the warning message
for k,v in warnings.items():
    warnings[k] = "[WARNING] "+v
# add stuff to rcParams
rcParams.update(warnings)

# add error-messages
errors = dict(
        E_failed_cons = "Failed to compute the consensus string.",
        E_zero_division = "Zero-division error encountered in '{0}' and '{1}'."
        )
for k,v in errors.items():
    errors[k] = "[ERROR] "+v
rcParams.update(errors)

# these are questions which need to be answered by the user
questions = dict(
        Q_errors_in_data = "There were errors in the input data. Do you want to exclude them?",
        )
for k,v in questions.items():
    questions[k] = "[?] "+v+" (y/n) "
rcParams.update(questions)

