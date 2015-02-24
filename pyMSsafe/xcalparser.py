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

This file handles the parsing of Xcalibur data structures.

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/cellzome/isobarquant
"""

import re
import copy
import numpy as np


# regex for the parsing of the FILTER field
REGEX = '(?P<analyser>[I|F]TMS) \+ (?P<centroid>[a-z]) [N|E]SI (det=(?P<det>\d*.\d*) )?'
REGEX += '(?P<xers>([a-zA-Z] )*)(?P<mode>Full|Z|SIM) ((?P<lock>lock) )?(?P<scan>ms[1-9]?) '
REGEX += '((?P<setmass1>\d*.\d*)@(?P<frag1>[a-z]*)(?P<energy1>\d*.\d*) )?'
REGEX += '((?P<setmass2>\d*.\d*)@(?P<frag2>[a-z]*)(?P<energy2>\d*.\d*) )?'
REGEX += '\[(?P<start>\d*.\d*)-(?P<end>\d*.\d*)\]'
RX_FILTER = re.compile(REGEX)

RX_REPEATS = re.compile('top (\d*) peaks')

RX_REJECT = re.compile('(\d*\.\d*)\s*(\d*\.\d*)\s*(\d*\.\d*)\s*(\d*\.\d*)\s*(\d*\.\d*)\s*(\d*\.\d*)')

rx_pump = re.compile('^(NC_Pump|LoadingPump)[\._](.*)')

endRX = '(?P<Start>\d*.\d*) *(?P<End>\d*.\d*)( *(?P<NCE>\d*.\d*))?( *(?P<CS>\d*))?( *(?P<Comment>.*))?'
rx_pol = '(?P<Polarity>Positive|Negative) *(?P<Mass>\d*.\d*) *' + endRX
rx_mz = '(?P<Mass>\d*.\d*) *(?P<Polarity>Positive|Negative) *' + endRX
rx_qex_incl_pol = re.compile(rx_pol)
rx_qex_incl_mz = re.compile(rx_mz)


def parseParameterList(data, wanted=''):
    """
    @brief Parses a parameter, no attempt is made to convert the value strings.
        The members are filtered by the wanted list.
    @param data <list>: list of parameter/value pairs
    @return a dictionary of the wanted values
    """
    # required variables
    params = {}
    sub = 0
    subname = ''
    subdic = {}
    subs = []

    # loop through the data
    for param_val in data:
        if param_val[0] == '' or param_val[0][:3] == '===':
            # signal that the next line will form a subdictionary
            sub = 1
            # copy an existing subdictionary to the appropriate name
            if subname:
                params[subname] = subdic.copy()

            if param_val[0][:3] == '===':
                # create new subdirectory
                subname = param_val[0][4:-6]
            else:
                subname = ''
            # clear the old data ready for the new subdictionary
            subdic = {}

        elif sub == 1 and subname == '':
            # this is the new subdictionary name
            subname = param_val[0]
            subs.append(subname)
        elif param_val[0][:-1] in wanted or not wanted:
            # parameter in the wanted list so append to the appropriate dictionary
            if subname:
                subdic[param_val[0][:-1]] = param_val[1]
            else:
                params[param_val[0][:-1]] = param_val[1]

    # add the last subdictionary
    if sub:
        params[subname] = subdic

    return params


def formatParamterList(unformatted):
    """
    @brief Parses instrument parameters, no attempt is made to convert the value strings.
    @param unformatted <string>: string of parameters and values
    @return a list of the split parameter - value pairs
    """

    # first split the data into lines and parameter:value pairs
    data = []
    for d in unformatted.splitlines():
        param_val = d.split(':')
        data.append(param_val)

    # process each paramter into tuples
    params = []
    for d in range(len(data)):
        if data[d][0] != '':
            # has parameter name
            par = data[d][0].strip()
            if len(data[d]) == 1:
                # no value given so check the parameter is not enabled etc
                for val in ['not enabled', 'enabled', 'disabled']:
                    pos = par.lower().find(val)
                    if pos != -1:
                        # found enabled
                        break
                if pos > 0:
                    # split the enabled etc from the parameter into the value
                    params.append((par[:pos], par[pos:]))
                else:
                    # leave as a parameter with no value
                    params.append((par, ''))
            else:
                # parameter and value added to the list
                params.append((par, data[d][1].strip()))
        else:
            # blank parameter only add if the next line should be a subdirectory (no value)
            if data[d + 1][0] == 'Segment 1 Information':
                params.append(('', ''))
            elif len(data[d + 1]) == 1:
                pass
            elif data[d + 1][1] != '':
                pass
            elif data[d + 1][0] == '(-) Global MS/MS Masses':
                pass
            else:
                params.append(('', ''))
    return params


def parseLTQ(data, skipEvents):
    """
    @brief Parses LTQ type instrument parameters, no attempt is made to convert the value strings.
    @param data <list>: list of parameter/value pairs
    @return a dictionary of the wanted values
    """
    params = {'main': {}}
    sub = 0
    subname = ''
    subs = []
    eventdic = {}
    eventlist = []
    special = 0
    reps = 0
    oldEventNum = 0
    globalParentMasses = []
    globalNonDataDependentSettings = {}
    globalRejectMasses = []
    microscans = []
    numInc = 0
    usincmasses = 0

    # do the formatting first
    paramlist = formatParamterList(data)

    # field names for the inclusion list
    inctitles = ['MS Mass', 'Start (min)', 'End (min)', 'MS FAIMS CV', 'MS Normalized Collision Energy',
                 'MS Charge State', 'MS Intensity Threshold', 'MS2 Mass', 'MS2 Normalized Collision Energy', 'Name']
    # lengths of the inclusion list fields
    inclens = [12, 15, 13, 11, 14, 10, 13, 12, 14, 8]

    msmstitles = ['MS Mass', 'Start (min)', 'End (min)', 'MS FAIMS CV',
                  'MS Normalized Collision Energy', 'Last Mass', 'Name']
    msmslens = [12, 15, 13, 11, 15, 10, 13]

    rejecttitles = ['Mass', 'Start (min)', 'End (min)', 'Mass', 'Start (min)', 'End (min)']

    event = ''
    order = ['main']
    for param_val in paramlist:
        if param_val[0] == '':
            # signal that the next line will form a subdictionary
            sub = 1
            # deal with the special subdictionaries
            if special == 1:
                params['repeats'] = reps
                if newEventNum > 1 and newEventNum not in skipEvents:
                    order.append(subname)
                    eventdic = updateeventdic(eventdic)
                    params[subname] = eventdic.copy()
                event = ''
            elif special == 2:
                order.append(subname)
                params[subname] = globalParentMasses
            elif special == 3:
                order.append(subname)
                params[subname] = globalNonDataDependentSettings
            elif special == 4:
                order.append(subname)
                params[subname] = globalRejectMasses
            elif special == 5:
                order.append(subname)
                params[subname] = microscans
            elif subname:
                order.append(subname)
                params[subname] = subdic.copy()
                if event in ('Chromatography Function', 'Dynamic Exclusion', 'Mass Tags'):
                    order.append(event)
                    params[event] = eventdic.copy()
                    event = ''

            # clear the old data ready for the new subdictionary
            subdic = {}
            subname = ''
            special = 0
        elif sub == 1 and subname == '':
            # this is the new subdictionary name
            subname = param_val[0]
            subs.append(subname)

            # some subdictionaries need special handling
            if subname == 'Scan Event Details':
                special = 1
            elif subname == 'Global Parent Masses':
                special = 2
                usincmasses = 0
                globalParentMasses = []
                numInc = 0
            elif subname == 'Global Non-Data Dependent Settings':
                special = 3
                globalNonDataDependentSettings = {}
            elif subname == 'Global Reject Masses':
                special = 4
                globalRejectMasses = []
                numInc = 0
            elif subname == 'Additional Microscans':
                special = 5
                microscans = []
                numInc = 0
            else:
                special = 0
        elif special == 1:
            # special: Scan Event Details
            try:
                # new event
                newEventNum = int(param_val[0])
                params['main']['num scan events'] = newEventNum
                params['num scan events'] = newEventNum
                if oldEventNum > 0 and oldEventNum not in skipEvents:
                    eventdic = updateeventdic(eventdic)
                    order.append(subname)
                    params[subname] = eventdic.copy()

                subname = 'Scan Event ' + param_val[0]
                event = newEventNum
                eventdic = dict(Type=param_val[1])
                oldEventNum = newEventNum
            except:
                if param_val[0][:10] == 'Scan Event':
                    # gives the repeats
                    rxm = RX_REPEATS.search(param_val[0])
                    reps = int(rxm.group(1))
                else:
                    eventdic[param_val[0].strip()] = param_val[1]
        elif special == 2:
            # special: Global Parent Masses
            numInc += 1
            try:
                if usincmasses:
                    globalParentMasses.append(formatGlobalMasses(param_val[0], inctitles, inclens))
                else:
                    masses = formatRejectMasses(param_val[0], rejecttitles)
                    for m in masses:
                        data = {'MS Mass': float(m['Mass']), 'Start (min)': m['Start (min)'],
                                'End (min)': m['End (min)'], 'MS FAIMS CV': '0', 'MS Normalized Collision Energy': '0',
                                'MS Charge State': '0', 'MS Intensity Threshold': '0', 'MS2 Mass': '0',
                                'MS2 Normalized Collision Energy': '0', 'Name': ''}
                        globalParentMasses.append(data)

            except ValueError:
                # finds which format of inclusion list is being used
                if param_val[0].find('FAIMS') > -1:
                    usincmasses += 1
        elif special == 3:
            # special: Global Non-Data Dependent Settings
            if param_val[0][-12:] == 'MS/MS Masses':
                # set-up the + or - MS/MS mass list
                parName = param_val[0]
                globalNonDataDependentSettings[parName] = []
            elif param_val[0] == '(none)':
                # deal with an empty set of masses
                globalNonDataDependentSettings[parName] = param_val[0]

            try:
                # try to read masses from the parameter data
                globalNonDataDependentSettings[parName].append(formatGlobalMasses(param_val[0], msmstitles, msmslens))
            except ValueError:
                # failure means that the data was the title line of the list
                pass
        elif special == 4:
            # special: Global Reject Masses
            numInc += 1
            try:
                # i = float(param_val[0][0:msmslens[0]])
                if not param_val[0].startswith('Mass'):
                    rejects = formatRejectMasses(param_val[0], rejecttitles)
                    globalRejectMasses += rejects
                pass
            except ValueError:
                pass
        elif special == 5:
            # special: Additional Microscans
            numInc += 1
            microscans.append(param_val[0].split())
            pass
        elif event == 'Additional Microscans':
            if param_val[0][:2] == 'MS':
                # accumulate microscan data to one list
                eventlist.append(param_val[0].strip())
            else:
                # end of microscan data so write list
                subdic[event] = eventlist
                event = ''
                # don't forget to add the next parameter
                subdic[param_val[0].strip()] = param_val[1].strip()
        elif event == 'Chromatography Function':
            if param_val[0] in ['Master scan event', 'Expected peak width', 'Min signal threshold',
                                'Method', 'Max area ratio to previous scan']:
                eventdic[param_val[0]] = param_val[1]
        elif event == 'Dynamic Exclusion':
            if param_val[0] in ['Repeat Count', 'Repeat Duration', 'Exclusion List Size',
                                'Exclusion Duration', 'Exclusion mass width relative to mass',
                                'Exclusion mass width relative to low (ppm)',
                                'Exclusion mass width relative to high (ppm)',
                                'Expiration']:
                eventdic[param_val[0]] = param_val[1]
                if param_val[0] == 'Expiration':
                    event = ''
        elif event == 'Mass Tags':
            if param_val[0] in ['Ratio Range (%)', 'Mass deltas']:
                eventdic[param_val[0]] = param_val[1]
        elif subname:
            if param_val[0] == 'Additional Microscans':
                # microscan data should be accumulated to a list
                event = 'Additional Microscans'
                eventlist = [param_val[1].strip()]
            elif param_val[0] == 'Neutral Loss Mass List':
                if param_val[1] == '(none)':
                    params['main'][param_val[0]] = param_val[1]
                else:
                    params['main'][param_val[0]] = param_val[1].strip()
            elif param_val[0] == 'Chromatography mode is ':
                event = 'Chromatography Function'
                eventdic = {'Chromatography mode': param_val[1]}
            elif param_val[0] == 'Dynamic exclusion ':
                event = 'Dynamic Exclusion'
                eventdic = {'Dynamic exclusion': param_val[1]}
            elif param_val[0] == 'Mass Tags data dependence ':
                event = 'Mass Tags'
                eventdic = {'Mass Tags data dependence': param_val[1]}
            else:
                # subdictionay being used so append to this
                subdic[param_val[0].strip()] = param_val[1].strip()

        else:
            if param_val[0] == 'Additional Microscans':
                # microscan data should be accumulated to a list
                event = 'Additional Microscans'
                eventlist = [param_val[1].strip()]
            else:
                # no subdictionary so append to the main dictionary
                params['main'][param_val[0].strip()] = param_val[1].strip()

    # add the last subdictionary
    if sub:
        order.append(subname)
        params[subname] = subdic

    return params, order


def string2int(value):
    try:
        val = int(value)
    except ValueError:
        val = -1
    return val


def string2float(value):
    try:
        val = float(value)
    except ValueError:
        val = -1
    return val


def parseQExative(data):
    """
    @brief Parses Q Exactive type instrument parameters, no attempt is made to convert the value strings.
    @param data <list>: list of parameter/value pairs
    @return a dictionary of the wanted values
    """
    special = 0
    subname = ''
    subname2 = ''
    subdic = {}
    params = {}
    order = []
    splits = []

    paramlist = formatExactiveList(data)

    for param_val in paramlist:
        if param_val[0] == 'General':
            continue

        if param_val[0] == '':
            # signal the end of a group and the beginning of a new group
            if subname and subdic:
                # has a subgroup and also parameters
                order.append(subname)
                params[subname] = copy.deepcopy(subdic)
            subname = ''
            subname2 = ''
            subdic = {}
            special = 0
        elif special == 1:
            if param_val[0] in ['related file', 'entries']:
                # file containing list data
                subdic[param_val[0]] = param_val[1]
                types = []
            elif types:
                try:
                    v = param_val[0]
                    s = splits[:]
                    float(v[:s[0]])
                    vals = []
                    t = 0
                    while len(s) > 0:
                        name = types[t][0]
                        if name in ['Start', 'Mass']:
                            vals.append(string2float(v[:s[0]]))
                        elif name in ['End']:
                            if v[:s[0]].strip() == 'end':
                                vals.append(string2float(runtime))
                            else:
                                vals.append(string2float(v[:s[0]]))
                        elif name in ['CS']:
                            vals.append(string2int(v[:s[0]]))
                        elif name in ['NCE']:
                            vals.append(string2int(v[:s[0] - 1]))
                        elif name in ['Comment']:
                            vals.append(v.strip())
                        else:
                            vals.append(v[:s[0]].strip())
                        trim = s.pop(0)
                        v = v[trim:]
                        t += 1

                    arraydata.put(pos, tuple(vals))
                    pos += 1

                    if pos == subdic['entries']:
                        # array filled with data so add to subdic
                        subdic['data'] = arraydata[:]
                        pass
                except ValueError:
                    # skip the units line
                    pass
            else:
                titles = re.findall('\s*\w*', param_val[0])
                splits = []
                for t in titles[:-1]:
                    splits.append(len(t))
                    ti = t.strip()
                    if ti in ['Polarity']:
                        types.append((ti, 'S10'))
                    elif ti in ['Comment']:
                        types.append((ti, 'S50'))
                        splits[-1] = 50
                    elif ti in ['NCE', 'CS']:
                        types.append((ti, 'int'))
                    elif ti in ['nCE', 'CS']:
                        types.append(('NCE', 'int'))
                    else:
                        types.append((ti, 'float'))

                arraydata = np.ndarray(subdic['entries'], dtype=types)
                pos = 0

        elif special == 2:
            if param_val[0] in ['entries']:
                # file containing list data
                subdic[param_val[0]] = param_val[1]
                subdic['data'] = []
            elif isinstance(param_val[0], dict):
                subdic['data'] = param_val[:]
                pass
            else:
                pass
        elif len(param_val) == 1:
            # name only parameter
            if subname:
                subname2 = param_val[0]
                subdic[subname2] = {}
            else:
                subname = param_val[0]
                if subname == 'LOCK MASSES':
                    special = 1
                elif subname == 'INCLUSION LIST':
                    special = 2
        elif subname2:
            # subset exists so append to this subset
            subdic[subname2][param_val[0]] = param_val[1]
            if param_val[0] == 'Method duration':
                runtime = param_val[1].split()[0]
                params['MS Run Time (min)'] = runtime
        elif subname:
            # subset exists so append to this subset
            subdic[param_val[0]] = param_val[1]
        else:
            params[param_val[0]] = param_val[1]

    # group scan methods into an experiment
    if 'EXPERIMENT' not in params:
        exptypes = {'FULL MS / DD-MS2 (TOPN)': ['Full MS', 'dd-MS2 / dd-SIM'],
                    'TARGETED-SIM / DD-MS2': ['SIM', 'dd-MS2'],
                    'FULL MA / AIF': [],
                    'FULL MS - SIM': [],
                    'TARGETED-SIM': ['SIM'],
                    'TARGETED-MS2': ['MS2']}
        params['EXPERIMENT'] = {}
        order.append('EXPERIMENT')
        for key in exptypes:
            if key in params:
                params['EXPERIMENT'][key] = params.pop(key)
                order.remove(key)
        pass

    # create run time info
    order.append('MS Run Time (min)')

    # create scan event data
    expts = params['EXPERIMENT'].keys()
    if len(expts) == 1:
        # single experiment so use the expected data
        exp = expts[0]
        se1 = formMSevent(params['EXPERIMENT'][exp][exptypes[exp][0]])
        se2 = formMSMSevent(params['EXPERIMENT'][exp][exptypes[exp][1]])
    elif expts == ['TARGETED-MS2', 'TARGETED-SIM'] or expts == ['TARGETED-SIM', 'TARGETED-MS2']:
        se1 = formMSevent(params['EXPERIMENT']['TARGETED-SIM']['SIM'])
        se2 = formMSMSevent(params['EXPERIMENT']['TARGETED-MS2']['MS2'])
    else:
        raise

    order.append('Scan Event 1')
    params['Scan Event 1'] = se1
    order.append('Scan Event 2')
    params['Scan Event 2'] = se2
    params['num scan events'] = 2
    order.append('units')
    params['units'] = {0: [], 1: []}
    params['units'][0].append({'highmz': se1['highmz'], 'lowmz': se1['lowmz'],
                               'resolution': se1['resolution'], 'scans': [1]})
    params['units'][1].append({'acttime': -1.0, 'activation': se2['activation'], 'resolution': se2['resolution'],
                               'energy': se2['energy'], 'isolation': se2['isolation'], 'lowmz': se2['lowmz'],
                               'colenergysteps': se2['ces'], 'colenergywidth': se2['cew'],
                               'scans': [x + 2 for x in range(se2['repeats'])]})
    params['activation'] = {1: ['HCD']}
    order.append('unit problems')
    params['unit problems'] = 0
    return params, order


def formatExactiveList(unformatted):
    """
    @brief Parses instrument parameters, no attempt is made to convert the value strings.
    @param unformatted <string>: string of parameters and values
    @return a list of the split parameter - value pairs
    """

    data = []
    inclist = []
    # remove unwanted special characters
    unformatted = unformatted.replace('\xb2', '2')
    unformatted = unformatted.replace('\x97', '-')
    # allow splitting for tune file
    lines = unformatted.replace('Tunefile', 'Tunefile ').split('\r\n')

    special = 0
    for l in lines:
        if l == '':
            # blank line: don't keep multiple blank lines
            if data[-1] != ['']:
                if special == 2 and inclist:
                    data.append(inclist)
                    inclist = []
                data.append([''])
            special = 0
        else:
            # parameter not blank line
            param_val = l.split('  ', 1)
            if special == 1:
                # treatment for Lock Mass data
                if l == '   (no file associated)':
                    data.append(['related file', 'none'])
                elif l == '   (no entries)':
                    data.append(['entries', 0])
                elif l.find('related') != -1:
                    param_val = l.strip().split(': ', 1)
                    data.append(param_val)
                elif l.find('entr') != -1:
                    tmp = l.strip().split(' ', 1)
                    data.append(['entries', int(tmp[0])])
                else:
                    data.append([l])
            elif special == 2:
                # deal with inclusion list data
                if l == '   (no file associated)':
                    data.append(['related file', 'none'])
                elif l == '   (no entries)':
                    data.append(['entries', 0])
                elif l.find('entr') != -1:
                    tmp = l.strip().split(' ', 1)
                    data.append(['entries', int(tmp[0])])
                elif l.strip()[:8] == 'Polarity':
                    # field names for the inclusion list
                    inctitles = ['Polarity', 'Mass', 'Start', 'End', 'NCE', 'CS', 'Comment']
                    # lengths of the inclusion list fields
                    inclens = [9, 10, 6, 7, 4, 4, 30]
                    inclist = []
                elif l.strip()[:4] == 'Mass':
                    # field names for the inclusion list
                    if l.find('Formula') != -1:
                        inctitles = ['Mass', 'Formula', 'Species', 'CS', 'Polarity', 'Start', 'End', 'NCE', 'Comment']
                        # lengths of the inclusion list fields
                        inclens = [10, 8, 8, 4, 9, 7, 7, 5, 30]
                    else:
                        inctitles = ['Mass', 'Polarity', 'Start', 'End', 'NCE', 'CS', 'Comment']
                        # lengths of the inclusion list fields
                        inclens = [10, 9, 6, 7, 4, 4, 30]
                    inclist = []
                elif l.find('tive') != -1:
                    # data for inclusion/exclusion lists
                    inclist.append(formatGlobalMasses(l, inctitles, inclens, 'CS'))
                else:
                    data.append([l])
                pass
            elif len(param_val) == 1:
                # only one entry this is title

                # check if title is for special data
                if l in ['ISOTOPIC PATTERN RECOGNITION', 'LOCK MASSES', 'NEUTRAL LOSSES', 'PRODUCT IONS', 'MASS TAGS']:
                    special = 1
                    data.append(param_val)
                elif l in ['EXCLUSION LIST', 'INCLUSION LIST']:
                    special = 2
                    data.append(param_val)
                else:
                    data.append([l.strip()])
                pass
            elif param_val[0] == '':
                # this is indented title
                data.append([l.strip()])
            else:
                # found param - val pair so use!
                param_val[1] = param_val[1].strip()

                data.append(param_val)

    # trim unwanted lines
    data.pop(0)
    data.pop(2)
    return data


def formMSMSevent(msms):
    """
    @brief takes Exactive "dd-MS2 / dd-SIM" parameters to make Scan Event 2
    @param msms <dictionary>: containing the MS scan parameters
    @return event <dictionary>: containing the LTQ type event parameters
    """

    event = {'group': 'MS2 event',
             'resolution': int(msms['Resolution'].replace(',', '')),
             'activation': 'HCD'}
    key = 'NCE'
    if 'NCE / stepped NCE' in msms:
        key = 'NCE / stepped NCE'

    try:
        event['energy'] = float(msms[key])
    except ValueError:
        event['energy'] = 0.0

    tmp = msms['Isolation window'].split()
    event['isolation'] = float(tmp[0])
    if msms['Fixed first mass'] == '?':
        event['lowmz'] = -1
    else:
        tmp = msms['Fixed first mass'].split()
        event['lowmz'] = float(tmp[0])
    if 'TopN' in msms:
        event['repeats'] = int(msms['TopN'])
    else:
        event['repeats'] = 1

    # to be updated when stepped colision energy data has been acquired
    event['ces'] = -1
    event['cew'] = -1
    return event


def formMSevent(fullms):
    """
    @brief takes Exactive "Full MS" parameters to make Scan Event 1
    @param fullms <dictionary>: containing the MS scan parameters
    @return event <dictionary>: containing the LTQ type event parameters
    """
    span = fullms['Scan range'].split()

    event = {'group': 'MS event',
             'resolution': int(fullms['Resolution'].replace(',', '')),
             'highmz': float(span[2]),
             'microscans': int(fullms['Microscans']), 'lowmz': float(span[0])}

    return event


def updateeventdic(eventdic):
    """
    @brief Updates values in the eventdic.
    @param eventdic <dictionary>: dictionary of current events
    @return eventdic <dictionary>: containing the updated values
    """
    ms2 = eventdic['Type'].find('MS/MS')
    ms3 = eventdic['Type'].find('MS3')
    ndd = eventdic['Type'].find('->oE')
    if ms2 > 0:
        # MS/MS event so extract group info
        loc2 = eventdic['Type'].find(' from')
        eventdic['group'] = eventdic['Type'][ms2:loc2]
    elif ms3 > 0:
        # MS3 event need to link it to the MS2 event it is derived from
        eventdic['group'] = eventdic['Type'][ms3:]
        loc = eventdic['Type'].find('(')
        loc2 = eventdic['Type'].find(')')
        scan = int(eventdic['Type'][loc + 1:loc2])

        # link the MS3 event to its MS2 event
        eventdic['ms2'] = scan
        eventdic['Activation Type'] += '3'
    elif ndd > 0:
        # non-data dependent scan
        eventdic['group'] = 'NDD event'
    elif 'Activation Type' in eventdic:
        # fragmentaion but not MS2 or MS3
        raise ValueError('Scan Event not compatible "%s"' % eventdic['Type'])
    else:
        eventdic['group'] = 'MS event'

    loc = eventdic['Type'].find('res=')
    if loc > 0:
        # found resolution setting
        loc2 = eventdic['Type'].find(' ', loc)
        eventdic['resolution'] = int(eventdic['Type'][loc + 4:loc2])

    # change the activation type if Multistage activation is enabled
    if 'Multistage activation is' in eventdic:
        eventdic['Activation Type'] += 'MSA'

    # add Wideband key if present in fragmentaion method
    if eventdic['Type'].find('Wideband') >= 0:
        eventdic['Wideband'] = 'enabled'
    return eventdic


def formatGlobalMasses(data, titles, lengths, charge=''):
    """
    @brief extracts the inclusion mass data from data string using the titles as
        dictionary keys and the lengths to split the fields
    @param data <string>: containing the inclusion mass data
    @param titles <list>: containing the dictionary keys
    @param lengths <list>: containing the number of characters in each field
    @return mass <dictionary>: containing the inclusion mass data
    """
    if charge:
        csKey = charge
    else:
        csKey = 'MS Charge State'
    start = 0
    mass = {}
    # data is in fixed width fields given by inclen
    for i in range(len(titles)):
        end = start + lengths[i]
        mass[titles[i]] = data[start:end].strip()
        if titles[i] == 'NCE' and mass[titles[i]][-1] == '%':
            mass[titles[i]] = mass[titles[i]][:-1].strip()
        start = end

    if csKey not in mass or mass[csKey] == '':
        mass[csKey] = '0'

    return mass


def formatRejectMasses(data, titles):  # , lengths):
    """
    @brief extracts the inclusion mass data from data string using the titles as dictionary keys
    @param data <string>: containing the inclusion mass data
    @param titles <list>: containing the dictionary keys
    @return masses <dictionary>: containing the inclusion mass data
    """
    masses = []
    mass = {}
    # data is in fixed width fields given by inclens
    grps = data.split()
    if 'none' in grps[0]:
        return masses

    for i in range(min(len(titles), len(grps))):
        if i % 3 == 0:
            if mass:
                masses.append(mass)
            mass = {}
        mass[titles[i]] = grps[i]
    masses.append(mass)
    return masses


def parseLCparameters(unformatted):
    """
    @brief Parses LC instrument parameters, no attempt is made to convert the value strings.
    @param unformatted <string>: string of parameters and values
    @return a list of the split parameter - value pairs
    """
    data = {}
    points = 0
    # split the lines and ittarate over the lines
    for l in unformatted.splitlines():
        # split paramter and value using '='
        param_val = l.split('=')
        par = param_val[0].strip()
        if len(param_val) == 2:
            # two values so add normally
            if par[-14:] == 'Profile_numpts':
                # special section for the pump flow conditions
                points = int(param_val[1])
                profile = par[:-7]
                prof = dict(numpts=points)
            elif par[3:] == 'Profile':
                # new layout with profile data on one line
                pts = param_val[1].split('|')
                prof = dict(pts=[])
                for p in pts:
                    try:
                        hplc = eval(p)
                        prof['pts'].append(dict(t=str(hplc[0]), t_unit='sec', flow=str(hplc[1]), flow_unit='nl/min'))
                    except:
                        continue
                prof['numpts'] = len(prof['pts'])
                data[par] = prof.copy()
            else:
                # normal parameter
                data[par] = param_val[1]
        elif points:
            # capture the pump data
            points -= 1
            val_unit = par[:-1].split('(')
            vals = val_unit[0].split(',')
            units = val_unit[1].split(',')
            # convert to dictionary
            dic = dict(t=vals[0], t_unit=units[0], flow=vals[1], flow_unit=units[1])
            try:
                prof['pts'].append(dic)
            except:
                prof['pts'] = [dic]

            # once all the points have been collected then add to the main dictionary
            if points == 0:
                data[profile] = prof.copy()

    return data


def parseDionexLCparameters(unformatted, pumps):
    """
    @brief Parses LC instrument parameters, no attempt is made to convert the value strings.
    @param unformatted <string>: string of parameters and values
    @return a list of the split parameter - value pairs
    """
    data = {'analytical': {'profile': []}, 'loading': {'profile': []}}
    time = -1
    # lines = unformatted.decode('latin-1').split('\n')
    lines = unformatted.split('\n')

    for li in lines:
        trim = li.strip()
        if not trim:
            continue
        if trim[0] == ';':
            # comment line: ignore
            # print 'comment\t' + trim
            pass
        elif trim[0] in '0123456789':
            # number so is a time point
            if time > -1:
                # new time point so save the old one
                for pump in profile:
                    if len(profile[pump]) > 1:
                        # data has been added to the pump profile so append to list
                        data[pump]['profile'].append(profile[pump].copy())

            time_param = trim.split(' ', 1)
            time = float(time_param[0])
            profile = {'analytical': {'time': time}, 'loading': {'time': time}}
            if len(time_param) > 1:
                param, val, pump = findParamValue(time_param[1].strip(), pumps)
                if val and pump:
                    if param == 'Flow':
                        flow_unit = val.split(' ', 1)
                        profile[pump]['flow'] = float(flow_unit[0])
                        profile[pump]['flow_unit'] = flow_unit[1]
                    elif param == '%B':
                        profile[pump]['b'] = float(val[:-1])
        else:
            # parameter line
            param, val, pump = findParamValue(trim, pumps)
            if val:
                if time > -1:
                    # data now the profile
                    if val and pump:
                        if param == 'Flow':
                            flow_unit = val.split(' ', 1)
                            profile[pump]['flow'] = float(flow_unit[0])
                            profile[pump]['flow_unit'] = flow_unit[1]
                        elif param == '%B':
                            profile[pump]['b'] = float(val[:-1])
                elif pump:
                    # this is a pump parameter in the main section so put in appropriate pump section
                    data[pump][param] = val
                else:
                    data[param] = val
    if time > -1 and len(profile) > 1:
        for pump in profile:
            data[pump]['points'] = len(data[pump]['profile'])
            if len(profile[pump]) > 1:
                # data has been added to the pump profile so append to list
                data[pump]['profile'].append(profile[pump].copy())

    return data


def findParamValue(data, pumps):
    """
    """
    val = ''
    pump = ''
    param_val = data.split('=', 1)
    param = param_val[0].strip()

    ispump = rx_pump.search(param)
    if ispump:
        gps = ispump.groups()
        pump = pumps[gps[0]]
        param = gps[1]

    if len(param_val) == 2:
        val = param_val[1].replace('[', '')
        val = val.replace(']', '')
        val = val.strip()
    return param, val, pump


def parseFilter(strFilter):
    """
    @brief Parses Filter string into its separate components using regular expression.
    @param strFilter <string>: string embedded values from Xcalibur
    @return a dictionary of the relevant values
    """
    rx_search = RX_FILTER.search(strFilter)
    if not rx_search:
        return {}
    return rx_search.groupdict()


def parseEvents(ms_param, skipEvents):
    """
    @brief Takes the Scan Events and creates the minimum set of Events to describe the MS method
    @param ms_param <dictionary>: containing details of the MS method
    @return events <dictionary>: containing the set of Events for the method.
    """
    scanevents = {}
    events = {}
    # collect the scan events into groups with the appropriate data
    for num in range(1, ms_param['num scan events'] + 1):
        if num in skipEvents:
            continue
        scanevent = ms_param['Scan Event %d' % num]
        if scanevent['group'] == 'MS event':
            if 'resolution' in scanevent:
                res = scanevent['resolution']
            else:
                res = -1
            data = {'highmz': scanevent['high mz'], 'lowmz': scanevent['low mz'], 'scans': [num], 'resolution': res}
            if 0 in scanevents:
                scanevents[0]['scans'].append(num)
            else:
                scanevents[0] = [data.copy()]
        elif scanevent['group'] == 'NDD event':
            params = scanevent['MS/MS'].split()
            data = {'activation': params[1], 'energy': float(params[3][:-1]), 'isolation': float(params[9]),
                    'scans': [num]}
            scanevents[0] = [data.copy()]
        else:
            if scanevent['Activation Type'] == 'HCD':
                lowmz = float(scanevent.get('FT first mass value', -1))
                if 'resolution' in scanevent:
                    res = scanevent['resolution']
                else:
                    res = -1
                # stepped col energy
                if ms_param['MS Detector Settings']['Stepped collision energy'] == 'enabled':
                    colsteps = int(ms_param['MS Detector Settings']['Number of collision energy steps'])
                    colwidth = float(ms_param['MS Detector Settings']['Collision energy width'])
                else:
                    colsteps = -1
                    colwidth = -1
            else:
                lowmz = -1
                res = -1
                colsteps = -1
                colwidth = -1

            data = {'activation': scanevent['Activation Type'], 'isolation': float(scanevent['Isolation Width']),
                    'energy': float(scanevent['Normalized Coll. Energy']), 'lowmz': lowmz, 'scans': [num],
                    'resolution': res, 'acttime': float(scanevent['Activation Time']), 'colenergysteps': colsteps,
                    'colenergywidth': colwidth}

            if 'ms2' in scanevent:
                # this is MS3 data and needs to be linked to the correct MS2 scan
                found = 0
                for unit in scanevents:
                    for scnevent in scanevents[unit]:
                        if scanevent['ms2'] in scnevent['scans']:
                            # found matching event
                            found = 1
                            break
                    if found:
                        break
            else:
                unit = events.setdefault(scanevent['group'], len(scanevents))

            scanevents.setdefault(unit, []).append(data.copy())

    # generate the minimum set of events and get activation list
    scanevents = colapseEvents(scanevents)
    acts = getactivation(scanevents)
    problems = 0

    if len(scanevents) > 2:
        # problem as too many events - only expect MS Event and one repeated MS/MS event
        problems = copy.deepcopy(scanevents)

        # first see if all activation methods the same
        base = 1

        for unit in range(2, len(scanevents)):
            errors = []
            for eventnum in range(len(scanevents[unit])):
                match, diffs = testeventmatch(scanevents[base][eventnum], scanevents[unit][eventnum])
                if not match:
                    errors.append(diffs)

            if errors:
                # find the event with the most spectra and treat as correct
                if len(scanevents[unit][eventnum]['scans']) > len(scanevents[base][eventnum]['scans']):
                    # use unit as best
                    best = unit
                    mrg = base
                else:
                    # use base as best
                    best = base
                    mrg = unit

                merged = scanevents[best][:]
                for eventnum in range(len(merged)):
                    merged[eventnum]['scans'] += scanevents[mrg][eventnum]['scans']
                    merged[eventnum]['scans'].sort()

                del scanevents[unit]
                del scanevents[base]
                scanevents[base] = merged[:]

        acts = getactivation(scanevents)

    # add extra scans for repeated events
    if ms_param['repeats']:
        # find the current max and the number of repeats already used
        maxscan = 0
        repsdone = 0
        for frag in scanevents[1]:
            maxscan = max(maxscan, max(frag['scans']))
            repsdone = len(frag['scans'])

        for rep in range(repsdone, ms_param['repeats']):
            for frag in scanevents[1]:
                maxscan += 1
                frag['scans'].append(maxscan)

    return scanevents, acts, problems


def testeventmatch(event1, event2):
    """
    @brief tests the correspondence between the 2 events
    @param event1 <dictionary>: containing the data for event 1
    @param event2 <dictionary>: containing the data for event 2
    @return <boolean>: truth of event matching
    @return diffs <dictionary>: highlighting the changed values
    """
    keys = ['activation', 'isolation', 'energy', 'lowmz', 'acttime']
    diffs = {}
    match = 0
    for k in keys:
        if event1[k] == event2[k]:
            match += 1
        else:
            diffs[k] = (event1[k], event2[k])

    return match == len(keys), diffs


def getactivation(scanevents):
    """
    @brief takes the raw set of Scan events and collapses them into the minimum set of events
    @param scanevents <dictionary>: containing all the event information
    """
    acts = {}
    for unit in range(1, len(scanevents)):
        for eventnum in scanevents[unit]:
            acts.setdefault(unit, []).append(eventnum['activation'])
    return acts


def colapseEvents(scanevents):
    """
    @brief takes the raw set of Scan events and collapses them into the minimum set of events
    @param scanevents <dictionary>: containing all the event information
    """
    okevents = [1]
    for unit in range(2, len(scanevents)):
        nomatch = 1
        lenev = len(scanevents[unit])
        for okunit in okevents:
            if lenev == len(scanevents[okunit]):
                matchfrag = 0
                for frag in range(lenev):
                    match, diffs = testeventmatch(scanevents[unit][frag], scanevents[okunit][frag])
                    if match:
                        matchfrag += 1
                if matchfrag == lenev:
                    # event matches so copy the scans into the appropriate methods
                    for frag in range(lenev):
                        scanevents[okunit][frag]['scans'] += scanevents[unit][frag]['scans']
                    # now delete from scanevents dictionary
                    del scanevents[unit]
                    nomatch = 0
                    break
        if nomatch:
            # new event so keep and add to list of okevents
            eventID = len(okevents) + 1
            okevents.append(eventID)
            if unit != eventID:
                scanevents[eventID] = scanevents[unit]
                del scanevents[unit]

    return scanevents
