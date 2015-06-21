"""
This module is part of the isobarQuant package,
written by Toby Mathieson and Gavain Sweetman
(c) 2015 Cellzome GmbH, a GSK Company, Meyerhofstrasse 1,
69117, Heidelberg, Germany.

The isobarQuant package processes data from
.raw files acquired on Thermo Scientific Orbitrap / QExactive
instrumentation working in  HCD / HCD or CID / HCD fragmentation modes.
It creates an .hdf5 file into which are later parsed the results from
Mascot searches. From these files protein groups are inferred and quantified.

This file mainly does exception handling.  Collecting all the exception classes
and formatting the output for reporting.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob/

"""

from pathlib import Path
import sys
import traceback


def trimTraceback():
    # gets the multiline format of the traceback and removes the last line
    tb = traceback.format_exc()
    end = tb[:-1].rfind('\n', )
    tb = tb[:end]

    return tb


def reformatException(error):
    exName = sys.exc_info()[0].__name__

    if exName == 'WindowsError':
        error.message = error.strerror

    elif len(error.args) > 1:
        # multiple args giving context
        error.message = error.args[0]
        addContext(error, error.args[1])

    return


def addContext(error, newContext):
    if not hasattr(error, 'context'):
        error.context = []
    error.context.append(newContext)


def exportError2File(error, file):
    """
    @brief takes the error object and writes the error data to the specified file
    @param error <Exception object>: containing the error data
    @param file <string>: containing the path for the error file
    @return:
    """

    filePath = Path(file)

    # create a dictionary of the error data
    errorDict = dict(code=1, error=error.message, type=sys.exc_info()[0].__name__, traceback=trimTraceback(),
                     repr=oneLineRepr(error))

    # add context only if present in the errror
    if hasattr(error, 'context'):
        errorDict['context'] = error.context

    # write to file
    fout = open(str(filePath), 'w')
    fout.write(str(errorDict))
    fout.close()


def multiLineRepr(error):
    """
    @brief turns the error object into multiple lines of text
    @param error <Exception object>: containing the error data
    @return outStr <string>: the formatted output
    """

    outStr = '%s' % trimTraceback()
    outStr += '\nError:    %s\nMessage:  %s' % (sys.exc_info()[0].__name__, error.message)
    if hasattr(error, 'context'):
        for c in error.context:
            outStr += '\nContext:  %s' % c

    return outStr


def oneLineRepr(error):
    """
    @brief turns the error object into multiple lines of text
    @param error <Exception object>: containing the error data
    @return outStr <string>: the formatted output
    """
    outStr = 'Error: %s, Message: %s' % (sys.exc_info()[0].__name__, error.message)
    if hasattr(error, 'context'):
        outStr += ', Context:  %s' % (' | '.join(error.context))

    return outStr


class czException(Exception):
    pass


class HDF5consistancyError(czException):
    pass


class HDFwritingError(czException):
    pass


class HDFmissingDataError(czException):
    pass


class HDFindexingError(czException):
    pass


class MGFprocessingError(czException):
    pass


class ImpFilesImportedError(czException):
    pass


class FileNotFoundException(czException):
    pass


class MissingParameter(czException):
    pass


class UsageError(czException):
    pass


class configError(czException):
    pass


class DBnoResults(czException):
    pass


class DBnoData(czException):
    pass


class DBextraData(czException):
    pass


class DBincompatibleDataType(czException):
    pass


class MascotModificationError(czException):
    pass


class MSdataConsistancyError(czException):
    pass


class XICprocessingError(czException):
    pass


class SpectraProcessingError(czException):
    pass


class FragmentMethodError(czException):
    pass


class QuantificationMethodError(czException):
    pass
