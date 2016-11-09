"""
This module is part of the isobarQuant package,
written by Toby Mathieson and Gavain Sweetman
(c) 2016 Cellzome GmbH, a GSK Company, Meyerhofstrasse 1,
69117, Heidelberg, Germany.

The isobarQuant package processes data from
.raw files acquired on Thermo Scientific Orbitrap / QExactive / Fusion
instrumentation working in  HCD / HCD or CID / HCD fragmentation modes.
It creates an .hdf5 file into which are later parsed the results from
Mascot searches. From these files protein groups are inferred and quantified.

This file mainly deals with the handling of the command line data and reading
configuration files converting parameter values as directed.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/protcode/isob

"""

from pathlib import Path
import ConfigParser
import types
import sys
import os
import re
from getopt import getopt
import ExceptionHandler as ExHa

regEx = re.compile('([A-Z_]*)')


class ConfigManager:
    """config container"""

    def __init__(self, configpath):
        self.parserpath = Path(configpath).absolute()
        config = ConfigParser.SafeConfigParser()
        config.readfp(open(str(self.parserpath)))

        # noinspection PyProtectedMember
        self.parameters = config._sections.copy()

        self.sectionConnector = '.'

        self.useDefaults = False

        self.immutableParamters = ['paramconversion', '__name__']
        self.validShortOptions = {'h': {'cmnd': 'helpRequested = True',
                                        'descr': 'request help text'},
                                  'd': {'cmnd': 'self.useDefaults = True',
                                        'descr': ('use default parameters (as given in %s)' % self.parserpath)}}

        self.mandatoryParamSection = 'runtime'
        self.adjustedOptionalParams = {}

        # read the environment variables
        self.environment = {}

        self.convertParameters()

    def get(self, section, parameter):
        return self.parameters[section][parameter]

    def scalePpmMda(self):

        for name in self.parameters:
            section = self.parameters[name]
            for param in section:
                if param in ['tolppm', 'largeppm', 'smallppm']:
                    section[param] /= 1000000
                elif param in ['tolmda', 'largemda', 'smallmda']:
                    section[param] /= 1000

    def convertParameters(self):
        converters = dict(eval=self.evalParameters, path=self.convert2path, int=self.convert2int,
                          list=self.convert2list, float=self.convert2float, bool=self.convert2boolian)
        for secName in self.parameters:

            section = self.parameters[secName]
            if 'paramconversion' in section:
                if isinstance(section['paramconversion'], basestring):
                    convDefs = eval(section['paramconversion'])
                else:
                    convDefs = section['paramconversion']

                section['paramconversion'] = convDefs

                for cd in convDefs:
                    if cd in converters:
                        converters[cd](secName, convDefs[cd])
                    else:
                        raise ExHa.configError('"%s" is not a valid definition for parameter conversion.' % cd)
            self.convertEnvironmentVariables(section)

    def convertEnvironmentVariables(self, section):
        # check for additional parameters with environment variables
        for param in section:
            if isinstance(section[param], basestring) and '$' in section[param]:
                section[param] = self.replaceEnvironmentVariables(section[param])

    def needsConverting(self, value):
        if isinstance(value, unicode) or isinstance(value, str):
            # is string and therefore needs converting
            if value.lower() == 'none':
                # value is none so convert to NoneType
                return 'convert2none'
            elif '$' in value:
                # this is an environment name
                return 'environment'
            else:
                # do full conversion
                return 'convert'
        # already converted
        return 'converted'

    def evalParameters(self, secName, conversionList):
        if len(conversionList) > 0:
            for param in conversionList:
                try:
                    action = self.needsConverting(self.parameters[secName][param])
                    if action == 'convert':
                        self.parameters[secName][param] = eval(self.parameters[secName][param])
                    elif action == 'convert2none':
                        self.parameters[secName][param] = None
                except:
                    raise ExHa.configError('"%s" in section "%s" can not be evaluated.' % (param, secName))

    def convert2list(self, secName, conversionList):
        # lists are assumed to be comma separated - in the command line with no spaces either
        if len(conversionList) > 0:
            for param in conversionList:
                if isinstance(self.parameters[secName][param], list):
                    # parameter already converted, don't repeat
                    continue
                try:
                    action = self.needsConverting(self.parameters[secName][param])
                    if action == 'convert':
                        if self.parameters[secName][param][0] == '[':
                            # bracketed list
                            self.parameters[secName][param] = eval(self.parameters[secName][param])
                        else:
                            # non-bracketed list: split on commas
                            tmp = self.parameters[secName][param].split(',')
                            for idx, t in enumerate(tmp):
                                try:
                                    tmp[idx] = eval(t)
                                except:
                                    tmp[idx] = t
                            self.parameters[secName][param] = tmp[:]

                    elif action == 'convert2none':
                        self.parameters[secName][param] = None
                except ExHa.czException, czEx:
                    ExHa.addContext(czEx, '"%s" in section "%s" can not be converted to list.' % (param, secName))
                    raise
                except:
                    raise ExHa.configError('"%s" in section "%s" can not be converted to list.' % (param, secName))

    def convert2int(self, secName, conversionList):
        if len(conversionList) > 0:
            for param in conversionList:
                try:
                    parameter = self.parameters[secName][param]
                    action = self.needsConverting(parameter)
                    if action == 'convert':
                        self.parameters[secName][param] = int(parameter)
                    elif action == 'convert2none':
                        self.parameters[secName][param] = None
                    elif action == 'environment':
                        self.parameters[secName][param] = self.replaceEnvironmentVariables(parameter, 'int')
                except:
                    raise ExHa.configError('"%s" in section "%s" can not be converted to int.' % (param, secName))

    def convert2float(self, secName, conversionList):
        if len(conversionList) > 0:
            for param in conversionList:
                try:
                    parameter = self.parameters[secName][param]
                    action = self.needsConverting(parameter)
                    if action == 'convert':
                        self.parameters[secName][param] = float(parameter)
                    elif action == 'convert2none':
                        self.parameters[secName][param] = None
                    elif action == 'environment':
                        self.parameters[secName][param] = self.replaceEnvironmentVariables(parameter, 'float')
                except:
                    raise ExHa.configError('"%s" in section "%s" can not be converted to float.' % (param, secName))

    def convert2path(self, secName, conversionList):
        if len(conversionList) > 0:
            for param in conversionList:
                try:
                    parameter = self.parameters[secName][param]
                    action = self.needsConverting(parameter)
                    if action == 'convert':
                        self.parameters[secName][param] = Path(parameter)
                    elif action == 'convert2none':
                        self.parameters[secName][param] = None
                    elif action == 'environment':
                        self.parameters[secName][param] = self.replaceEnvironmentVariables(parameter, 'path')
                except:
                    raise ExHa.configError('"%s" in section "%s" can not be converted to Path.' % (param, secName))

    def convert2boolian(self, secName, conversionList):
        if len(conversionList) > 0:
            for param in conversionList:
                try:
                    action = self.needsConverting(self.parameters[secName][param])
                    if action == 'convert':
                        value = self.parameters[secName][param].lower()
                        if value in ['yes', 'true', '1']:
                            self.parameters[secName][param] = True
                        elif value in ['no', 'false', '0']:
                            self.parameters[secName][param] = False
                        else:
                            raise ExHa.configError('')
                    elif action == 'convert2none':
                        self.parameters[secName][param] = None
                    elif action == 'environment':
                        self.parameters[secName][param] = self.replaceEnvironmentVariables(
                            self.parameters[secName][param], 'bool')
                except:
                    raise ExHa.configError('"%s" in section "%s" can not be converted to Boolian.' % (param, secName))

    def replaceEnvironmentVariables(self, value, convertTo=0):
        environ = self.environment

        data = value.split('$')
        # if starts with $ get an empty string as first element

        for i in range(1, len(data)):
            name = regEx.search(data[i]).groups()[0]

            if name in environ:
                data[i] = data[i].replace(name, str(environ[name]))
            else:
                raise ExHa.configError('"%s" not found in environment variables' % name)

        x = ''.join(data)

        if convertTo == 'path':
            rtn = Path(''.join(data))
        elif convertTo == 'int':
            rtn = int(''.join(data))
        elif convertTo == 'float':
            rtn = float(''.join(data))
        elif convertTo == 'eval':
            rtn = eval(''.join(data))
        else:
            rtn = ''.join(data)

        return rtn

    def convertConfig(self):
        """
        @brief converts the config parmeters dictionary to a list of dictionaries with keys set, parameter & value
        """

        configList = []

        for set in self.parameters:
            if set in ['production', 'test', 'dev']:
                # skip the mode sets as these contain passwords
                continue
            setdata = self.parameters[set]
            for parameter in setdata:
                if parameter in ('__name__', 'paramconversion'):
                    # skip the name parameter as this doesn't contain data
                    continue
                configList.append(dict(set=set, parameter=parameter, value=str(setdata[parameter])))

        return configList

    def evaluateCommandLineArgs(self, clArgs):

        helpRequested = False
        optionalParamsAdjusted = False

        callerPathTokens = clArgs[0].split('/')
        calledFrom = callerPathTokens[-1]

        usageTxt = self.contructUsageTxt(calledFrom)

        usingDefaultsTxt = 'DEFAULT-FLAG HAS BEEN SET -> ALL OTHER COMMAND LINE ARGUMENTS WILL BE IGNORED' \
                           '\n\tUSING DEFAULT PARAMETERS:\n'

        mandatoryParamsCheckList = {}
        mandatoryParamsTxt = 'MANDATORY PARAMETERS:\n'
        adjustedOptionalParamsTxt = 'ADJUSTED OPTIONAL PARAMETERS:\n'
        parametersAcceptedTxt = 'RUNNING %s WITH THE FOLLOWING PARAMETERS:'
        missingMandatoryParamsTxt = 'NOT ALL MANDATORY PARAMETERS HAVE BEEN PROVIDED.'

        mandatoryParams = self.parameters[self.mandatoryParamSection]

        helpTxt = self.contructHelpTxt(mandatoryParams)

        if mandatoryParams:
            baseString = '\t%-' + str(self.mandatoryParamLength + 2) + 's %s\n'
            for mp in mandatoryParams:
                if mp not in self.immutableParamters:
                    usingDefaultsTxt += baseString % (mp + ':', str(mandatoryParams[mp]))
                    mandatoryParamsCheckList[mp] = 'NOT PROVIDED'

        # extract section names from parameters object
        paramSecNames = self.parameters.keys()

        # prepare a list of all parameters in the cfg file as a list of accepted parameters from the command line
        try:
            longOptsMaxLength = 0
            longOpts = []
            for secName in paramSecNames:
                currSec = self.parameters[secName]
                for paramName in currSec:
                    # some parameters can not be configured via the command line
                    if paramName not in self.immutableParamters:
                        loName = '%s%s%s=' % (secName, self.sectionConnector, paramName)
                        if len(loName) > longOptsMaxLength:
                            longOptsMaxLength = len(loName)
                        longOpts.append(loName)
                        #   for convenience, mandatory parameters can be provided without specifying the section
                        if secName == self.mandatoryParamSection:
                            longOpts.append('%s=' % paramName)
            self.longOptsMaxLength = longOptsMaxLength

        except Exception as genEx:
            ExHa.addContext(genEx, 'Error in preparing the options list from the config parameters')
            raise

        shortOpts = []

        for vso in self.validShortOptions:
            longOpts.append('%s' % vso)

        opts, args = getopt(clArgs[1:], shortOpts, longOpts)
        printLines = []

        if len(opts) > 0:
            # scan short option first to check whether default values should be used or help is requested
            for opt in opts:
                optName = opt[0]
                if len(optName) == 3 and optName[-1] in self.validShortOptions:
                    exec self.validShortOptions[optName[-1]]['cmnd']

            if helpRequested:
                # if help flag has been set, display usage and help information
                print usageTxt
                print helpTxt
                sys.exit()
            else:
                if self.useDefaults:
                    # if default flag set, retrieve all parameters from the cfg file and ignore other command line
                    # arguments report the status of the mandatory parameters and return correct state
                    printLines = [usingDefaultsTxt]
                else:
                    mandatoryString = '\t%-' + str(self.mandatoryParamLength + 2) + 's %s\n'
                    adjustedString = '\t%-' + str(self.longOptsMaxLength + 2) + 's %s\n'
                    for opt in opts:
                        optName = opt[0]
                        optValue = opt[1]
                        if len(optName) > 2:
                            if self.sectionConnector in optName:
                                optNameTokens = optName.split(self.sectionConnector)
                                secName = optNameTokens[0][2:]
                                paramName = optNameTokens[1]
                                self.parameters[secName][paramName] = optValue
                            else:
                                secName = self.mandatoryParamSection
                                paramName = optName[2:]
                                self.parameters[secName][paramName] = optValue

                            #   if the current option refers to a mandatory parameter, check this parameter as being set
                            if secName == self.mandatoryParamSection and optValue:
                                mandatoryParamsCheckList[paramName] = optValue
                            else:
                                #   otherwise report that an optional parameter has been adjusted
                                adjustedOptionalParamsTxt += adjustedString % (optName + ':', str(optValue))
                                optionalParamsAdjusted = True
                                if secName in self.adjustedOptionalParams:
                                    self.adjustedOptionalParams[secName][paramName] = optValue
                                else:
                                    self.adjustedOptionalParams[secName] = {paramName: optValue}

                    #   create the output text that reports the status of the mandatory parameters
                    if mandatoryParams:
                        for mp in mandatoryParamsCheckList:
                            if mp not in self.immutableParamters:
                                mandatoryParamsTxt += mandatoryString % (mp + ':', str(mandatoryParamsCheckList[mp]))

                    #   if not all mandatory parameters have been supplied and default-flag has not been set
                    #   report which mandatory parameters have been set and which haven't and return error state
                    if 'NOT PROVIDED' in mandatoryParamsCheckList.values():
                        raise ExHa.MissingParameter('Error not all manditatory parameters were given: ' +
                                                    mandatoryParamsTxt)
                    else:
                        #   if default-flag has not been set and all mandatory parameters HAVE been supplied
                        #   report the status of the mandatory and the adjusted optional parameters
                        printLines = [parametersAcceptedTxt % sys.argv[0],  mandatoryParamsTxt]
                        if optionalParamsAdjusted:
                            printLines.append(adjustedOptionalParamsTxt)

        # if no parameters have been supplied at all, print usage information
        else:
            raise ExHa.UsageError('No parameters were given:', usageTxt)

        # finally, convert adjusted entries in the config object to their corresponding data type
        #   (via the command line, they are retrieved as strings) and return correct state
        self.convertParameters()
        try:
            screenLevel = self.parameters['logging']['screenlevel'].lower()
        except KeyError:
            screenLevel = 'none'

        if screenLevel in ['none', 'info', 'debug']:
            for line in printLines:
                print line
            sys.stdout.flush()

        return args

    def evaluateCommandLineArgsOLD(self, clArgs):

        helpRequested = False
        optionalParamsAdjusted = False

        callerPathTokens = clArgs[0].split('/')
        calledFrom = callerPathTokens[-1]

        usageTxt = self.contructUsageTxt(calledFrom)

        usingDefaultsTxt = 'DEFAULT-FLAG HAS BEEN SET -> ALL OTHER COMMAND LINE ARGUMENTS WILL BE IGNORED' \
                           '\n\tUSING DEFAULT PARAMETERS:\n'

        mandatoryParamsCheckList = {}
        mandatoryParamsTxt = 'MANDATORY PARAMETERS:\n'
        adjustedOptionalParamsTxt = 'ADJUSTED OPTIONAL PARAMETERS:\n'
        parametersAcceptedTxt = 'RUNNING %s WITH THE FOLLOWING PARAMETERS:'
        missingMandatoryParamsTxt = 'NOT ALL MANDATORY PARAMETERS HAVE BEEN PROVIDED.'

        mandatoryParams = self.parameters[self.mandatoryParamSection]

        helpTxt = self.contructHelpTxt(mandatoryParams)

        if mandatoryParams:
            baseString = '\t%-' + str(self.mandatoryParamLength + 2) + 's %s\n'
            for mp in mandatoryParams:
                if mp not in self.immutableParamters:
                    usingDefaultsTxt += baseString % (mp + ':', str(mandatoryParams[mp]))
                    mandatoryParamsCheckList[mp] = 'NOT PROVIDED'

        # extract section names from parameters object
        paramSecNames = self.parameters.keys()

        # prepare a list of all parameters in the cfg file as a list of accepted parameters from the command line
        try:
            longOptsMaxLength = 0
            longOpts = []
            for secName in paramSecNames:
                currSec = self.parameters[secName]
                for paramName in currSec:
                    # some parameters can not be configured via the command line
                    if paramName not in self.immutableParamters:
                        loName = '%s%s%s=' % (secName, self.sectionConnector, paramName)
                        if len(loName) > longOptsMaxLength:
                            longOptsMaxLength = len(loName)
                        longOpts.append(loName)
                        #   for convenience, mandatory parameters can be provided without specifying the section
                        if secName == self.mandatoryParamSection:
                            longOpts.append('%s=' % paramName)
            self.longOptsMaxLength = longOptsMaxLength

        except Exception as genEx:
            ExHa.addContext(genEx, 'Error in preparing the options list from the config parameters')
            raise

        shortOpts = []

        for vso in self.validShortOptions:
            longOpts.append('%s' % vso)

        opts, args = getopt(clArgs[1:], shortOpts, longOpts)

        if len(opts) > 0:
            # scan short option first to check whether default values should be used or help is requested
            for opt in opts:
                optName = opt[0]
                if len(optName) == 3 and optName[-1] in self.validShortOptions:
                    exec self.validShortOptions[optName[-1]]['cmnd']

            if not helpRequested:
                if not self.useDefaults:
                    mandatoryString = '\t%-' + str(self.mandatoryParamLength + 2) + 's %s\n'
                    adjustedString = '\t%-' + str(self.longOptsMaxLength + 2) + 's %s\n'
                    for opt in opts:
                        optName = opt[0]
                        optValue = opt[1]
                        if len(optName) > 2:
                            if self.sectionConnector in optName:
                                optNameTokens = optName.split(self.sectionConnector)
                                secName = optNameTokens[0][2:]
                                paramName = optNameTokens[1]
                                self.parameters[secName][paramName] = optValue
                            else:
                                secName = self.mandatoryParamSection
                                paramName = optName[2:]
                                self.parameters[secName][paramName] = optValue

                            #   if the current option refers to a mandatory parameter, check this parameter as being set
                            if secName == self.mandatoryParamSection and optValue:
                                mandatoryParamsCheckList[paramName] = optValue
                            else:
                                #   otherwise report that an optional parameter has been adjusted
                                adjustedOptionalParamsTxt += adjustedString % (optName + ':', str(optValue))
                                optionalParamsAdjusted = True
                                if secName in self.adjustedOptionalParams:
                                    self.adjustedOptionalParams[secName][paramName] = optValue
                                else:
                                    self.adjustedOptionalParams[secName] = {paramName: optValue}

                    #   create the output text that reports the status of the mandatory parameters
                    if mandatoryParams:
                        for mp in mandatoryParamsCheckList:
                            if mp not in self.immutableParamters:
                                mandatoryParamsTxt += mandatoryString % (mp + ':', str(mandatoryParamsCheckList[mp]))

                    #   if not all mandatory parameters have been supplied and default-flag has not been set
                    #   report which mandatory parameters have been set and which haven't and return error state
                    if 'NOT PROVIDED' in mandatoryParamsCheckList.values():
                        raise ExHa.MissingParameter(
                            'Error not all manditatory parameters were given: ' + mandatoryParamsTxt)
                    else:
                        #   if default-flag has not been set and all mandatory parameters HAVE been supplied
                        #   report the status of the mandatory and the adjusted optional parameters
                        print(parametersAcceptedTxt % sys.argv[0])
                        print mandatoryParamsTxt
                        if optionalParamsAdjusted:
                            print adjustedOptionalParamsTxt
                        sys.stdout.flush()

                # if default flag set, retrieve all parameters from the cfg file and ignore other command line arguments
                # report the status of the mandatory parameters and return correct state
                else:
                    print usingDefaultsTxt
                    sys.stdout.flush()

            #   if help flag has been set, display usage and help information
            else:
                print(usageTxt)
                print(helpTxt)
                sys.exit()

        # if no parameters have been supplied at all, print usage information
        else:
            raise ExHa.UsageError('No parameters were given:', usageTxt)

        # finally, convert adjusted entries in the config object to their corresponding data type
        #   (via the command line, they are retrieved as strings) and return correct state
        self.convertParameters()

        return args

    def contructUsageTxt(self, calledFrom):
        usageTxt = 'USAGE:\n\t%s [any parameter given in %s in the format \n\t\t\t--<section name>'
        usageTxt += '%s<parameter name> <parameter value>]'
        usageTxt = usageTxt % (calledFrom, self.parserpath, self.sectionConnector)
        for shortOpt in self.validShortOptions:
            usageTxt += '\n\t%s --%s\t%s' % (calledFrom, shortOpt, self.validShortOptions[shortOpt]['descr'])

        return usageTxt

    def contructHelpTxt(self, mandatoryParams):
        helpTxt = '\nHELP:'
        helpTxt += '\n\tAny parameter given in %s can be adjusted with the following syntax:' % self.parserpath
        helpTxt += '\n\t\t--<section name>%s<parameter name> <parameter value>' % self.sectionConnector

        if mandatoryParams:
            helpTxt += '\n\n\tIf the default flag is not set, at least all mandatory parameters have to be provided.'
            helpTxt += '\n\tThese are:'
            maxLength = 0
            for mp in mandatoryParams:
                if mp not in self.immutableParamters:
                    helpTxt += '\n\t\t%s' % mp
                    if len(mp) > maxLength:
                        maxLength = len(mp)

            helpTxt += '\n\n\tIf you set the default flag all parameters will be used as given in %s. ' \
                       '\n\tThe mandatory parameters will be set to:' % self.parserpath

            baseString = '\n\t\t%-' + str(maxLength + 2) + 's %s'
            for mp in mandatoryParams:
                if mp not in self.immutableParamters:
                    helpTxt += baseString % (mp + ':', mandatoryParams[mp])
            self.mandatoryParamLength = maxLength
        return helpTxt


class pyMSsafeConfigManager(ConfigManager):
    def formatFragMeths(self):

        tmp = self.parameters['general']['fragmeths']
        allmeths = []
        meths = []
        for t in tmp:
            meth = []
            meths = t.split(',')
            for m in meths:
                act, use = m.split('-')
                meth.append((act, use))
            allmeths.append(meth[:])

        return allmeths
