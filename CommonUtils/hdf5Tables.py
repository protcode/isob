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

This file holds all the table definitions used in hdf5 files for the MS.hdf5
and results.hdf5 files.  It also contains the descriptions of the table
columns that are copied into the tables as [table object].attrs.columnDesc

isobarQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

A copy of the license should have been part of the
download. Alternatively it can be obtained here :
https://github.com/cellzome/isobarquant
"""

import tables as tables


class Spectra(tables.IsDescription):
    """
    @brief HDF5 table definition for Spectra table: holds spectrum types and RTs
    """
    spec_id = tables.Int32Col(pos=0)
    rt = tables.Float32Col(pos=1)
    type = tables.StringCol(5, pos=2)

Spectra_attr = ['COL spec_id: Spectrum identifier, corresponds to the Xcalibur spectrum number.',
                'COL rt: The retention time of the spectrum.',
                'COL type: Type of spectrum MS1, MS2 or MS3.']


class MSMSheader(tables.IsDescription):
    """
    @brief HDF5 table definition for MSMSjeader table: Holds all data for MS/MS events
    """
    spec_id = tables.Int32Col(pos=0)
    unit_id = tables.Int32Col(pos=1)
    order = tables.Int32Col(pos=2)
    scanevent = tables.Int32Col(pos=3)
    precmz = tables.Float32Col(pos=4)
    precmzraw = tables.Float32Col(pos=5)
    precinten = tables.Float32Col(pos=6)
    peaktype = tables.StringCol(10, pos=7)
    monomz = tables.Float32Col(pos=8)
    precmz_xic = tables.Float32Col(pos=9)
    precmz_surv = tables.Float32Col(pos=10)
    rt = tables.Float32Col(pos=11)
    rtapex = tables.Float32Col(pos=12)
    scanapex = tables.Int32Col(pos=13)
    charge = tables.Int32Col(pos=14)
    precsd = tables.Float32Col(pos=15)
    id_spec = tables.Int32Col(pos=16)
    quan_spec = tables.Int32Col(pos=17)
    survey_spec = tables.Int32Col(pos=18)
    inten = tables.Float32Col(pos=19)
    area = tables.Float32Col(pos=20)
    setmass = tables.Float32Col(pos=21)
    fwhm = tables.Float32Col(pos=22)
    c13iso = tables.Int32Col(pos=23)
    s2i = tables.Float32Col(pos=24)
    s2iearly = tables.Float32Col(pos=25)
    s2ilate = tables.Float32Col(pos=26)
    c12 = tables.Float32Col(pos=27)
    c12early = tables.Float32Col(pos=28)
    c12late = tables.Float32Col(pos=29)
    thresh = tables.Float32Col(pos=30)
    threshearly = tables.Float32Col(pos=31)
    threshlate = tables.Float32Col(pos=32)
    fragmeth = tables.StringCol(5, pos=33)
    fragenergy = tables.Float32Col(pos=34)
    sumrepint = tables.Float32Col(pos=35)
    sumreparea = tables.Float32Col(pos=36)
    smooth = tables.Float32Col(pos=37)
    integ = tables.Float32Col(pos=38)

MSMSheader_attr = ['COL spec_id: Spectrum identifier; corresponds to the Xcalibur spectrum number.',
                   'COL unit_id: Unit identifier: The unit groups spectra from the same precursor so that they can be processed together.',
                   'COL order: Order in which the spectra in a unit were acquired.',
                   'COL scanevent: The event number assigned by Xcalibur.',
                   'COL precmz: The m/z for the precursor ion using XIC derived data if possible.  If recalibration has been performed this is the recalibrated m/z',
                   'COL precmzraw: The m/z for the precursor ion using XIC derived data if possible.',
                   'COL precinten: The intensity value for the precursor ion. This will be the apex intensity of the XIC data where possible or the interpolated intensity of the precursor ion at the time of the MS/MS event.',
                   'COL peaktype: The source of the precursor ion data. Is One of XICclust - XIC data with full cluster matching, XICsingle - XIC data with only 12C ion matching, survey - interpolated data from MS1 spectra.',
                   'COL monomz: The Xcalibur monoisotopicMZ value - a m/z value for the precursor ion.',
                   'COL precmz_xic: The m/z for the precursor ion derived from XIC data.  Intensity weighted mean of the XIC data points with intensity > half the maximum intensity.',
                   'COL precmz_surv:  The m/z for the precursor ion derived from MS1 spectra data.',
                   'COL rt: The retention time of the spectrum.',
                   'COL rtapex: The retention time of the matching XIC peak apex.',
                   'COL scanapex: The spectrum number for the XIC peak apex.',
                   'COL charge: The charge state of the precursor ion.',
                   'COL precsd: The standard deviation of the data point used to calculate the precmz_xic.',
                   'COL id_spec: The spec_id of the spectrum where the sequence data is located.  When units have multiple spectra one can contain the sequence information and the other the quantification data.',
                   'COL quan_spec: The spec_id of the spectrum where the quantification data is located.  When units have multiple spectra one can contain the sequence information and the other the quantification data.',
                   'COL survey_spec: Spectrum identifier for the MS1 spectrum prior to the MS/MS event.',
                   'COL inten: The apex ion intensity of any linked XIC peak.',
                   'COL area: The maximum ion area of any linked XIC peak.',
                   'COL setmass: The m/z set by Xcalibur as the center of the isolation window.',
                   'COL fwhm: The width of the XIC peak at half maximum intensity.',
                   'COL c13iso: The number of 13C isotopes identified in the isolation window.',
                   'COL s2i: The S2I value for the MS/MS event, this is interpolated from the two flanking MS1 spectra.',
                   'COL s2iearly: The S2I value for the precursor ion measured in the MS1 event prior to the MS/MS event.',
                   'COL s2ilate: The S2I value for the precursor ion measured in the MS1 event after to the MS/MS event.',
                   'COL c12: The intensity of the precursor 12C ion, this is interpolated from the two flanking MS1 spectra.',
                   'COL c12early: The 12C intensity value for the precursor ion measured in the MS1 event prior to the MS/MS event.',
                   'COL c12late: The 12C intensity value for the precursor ion measured in the MS1 event after to the MS/MS event.',
                   'COL thresh: The Xcalibur intensity cut off (best estimation of noise) level for the MS/MS event, this is interpolated from the two flanking MS1 spectra.',
                   'COL threshearly: The Xcalibur intensity cut off for the precursor ion measured in the MS1 event prior to the MS/MS event.',
                   'COL threshlate: The Xcalibur intensity cut off for the precursor ion measured in the MS1 event after to the MS/MS event.',
                   'COL fragmeth: The fragmentation method used in this MS/MS event.',
                   'COL fragenergy: The energy used for fragmentation for this MS/MS event.',
                   'COL sumrepint: Sum of the intensities of all identified reporter ions.',
                   'COL sumreparea: Sum of the areas of all identified reporter ions.',
                   'COL smooth: The maximum smoothed intensity for the XIC peak.',
                   'COL integ: The integrated ion intensity for ions with intensities > half the maximum intensity the XIC peak.']


class SpecParams(tables.IsDescription):
    """
    @brief HDF5 table definition for Spectrum Parameters table:  Holds Xcalibur parameter data
    """
    spec_id = tables.Int32Col(pos=0)
    set = tables.StringCol(30, pos=1)
    parameter = tables.StringCol(60, pos=2)
    value = tables.StringCol(200, pos=3)

SpecParams_attr = ['COL spec_id: Spectrum identifier; corresponds to the Xcalibur spectrum number.',
                   'COL set: The Xcalibur group of parameters: either "filter" or "Extra".',
                   'COL parameter: The parameter name.',
                   'COL value: The parameter value.']


class XICbins(tables.IsDescription):
    """
    @brief HDF5 table definition for raw XIC ion data:  Holds MS1 data
    """
    bin = tables.Int32Col(pos=0)
    specid = tables.Int32Col(pos=1)
    rt = tables.Float32Col(pos=2)
    mz = tables.Float32Col(pos=3)
    inten = tables.Float32Col(pos=4)
    area = tables.Float32Col(pos=5)
    mzraw = tables.Float32Col(pos=6)

XICbins_attr = ['COL bin: The m/z bin, this is calculated from the m/z using the binsPerDa parameter.',
                'COL specid: Spectrum identifier; corresponds to the Xcalibur spectrum number.',
                'COL rt: The retention time of the spectrum.',
                'COL mz: The m/z of the ion, after recalibration if performed.',
                'COL inten: The intensity of the ion.',
                'COL area: The ion area, calculated by integrating the intensities of the ion.',
                'COL mzraw: The m/z of the ion.']


class Parameters(tables.IsDescription):
    """
    @brief HDF5 table definition for Parameters table
    """
    set = tables.StringCol(30, pos=0)  # from the call used to get the parameter
    subset = tables.StringCol(30, pos=1)  # subset information '' or name from Xcalibur
    parameter = tables.StringCol(60, pos=2)  # parameter name
    value = tables.StringCol(100, pos=3)  # parameter value as a string

Parameters_attr = ['COL set: What part of the system the parameters are from, e.g. LC, Instrument, Tune.',
                   'COL subset: The grouping of the parameters e.g. analytical, EXPERIMENT, Scan Event.',
                   'COL parameter: The name of the parameter.',
                   'COL value: The value of the parameter.']


class Ions(tables.IsDescription):
    """
    @brief HDF5 table definition for Ions table. Holds MS/MS ion data
    """
    spec_id = tables.Int32Col(pos=0)
    mz = tables.Float32Col(pos=1)
    inten = tables.Float32Col(pos=2)
    mzraw = tables.Float32Col(pos=3)
    rt = tables.Float32Col(pos=4)

Ions_attr = ['COL spec_id: Spectrum identifier; corresponds to the Xcalibur spectrum number.',
             'COL mz: The m/z for the ion, after recalibration if performed.',
             'COL inten: The intensity of the MS/MS ion.',
             'COL mzraw: The m/z for the ion.',
             'COL rt: The retention time of the spectrum.']


class DeconvIons(tables.IsDescription):
    """
    @brief HDF5 table definition for DeconvIons table:  Holds deconvoluted MS/MS spectra data
    """
    spec_id = tables.Int32Col(pos=0)
    mz = tables.Float32Col(pos=1)
    inten = tables.Float32Col(pos=2)
    charge = tables.StringCol(10, pos=3)

DeconvIons_attr = ['COL spec_id: Spectrum identifier; corresponds to the Xcalibur spectrum number.',
                   'COL mz: The m/z of the MS/MS ion.',
                   'COL inten: The intensity of the MS/MS ion.',
                   'COL charge: The charge state of the MS/MS ion.']


class Profile(tables.IsDescription):
    """
    @brief HDF5 table definition for LC gradient data table
    """
    id = tables.StringCol(6, pos=0)
    time = tables.Int32Col(pos=1)
    time_unit = tables.StringCol(3, pos=2)
    flow = tables.Float32Col(pos=3)
    flow_unit = tables.StringCol(6, pos=4)

Profile_attr = ['COL id: The name of the LC and solvent container.',
                'COL time: The time point.',
                'COL time_unit: The unit of time used in time.',
                'COL flow: The solvent flow rate.',
                'COL flow_unit: The unit of flow used in flow']


class Dionex(tables.IsDescription):
    """
    @brief HDF5 table definition for LC gradient data table for Dionex LC
    """
    id = tables.StringCol(12, pos=0)
    time = tables.Float32Col(pos=1)
    time_unit = tables.StringCol(3, pos=2)
    percent_b = tables.Float32Col(pos=3)
    flow = tables.Float32Col(pos=4)
    flow_unit = tables.StringCol(6, pos=5)

Dionex_attr = ['COL id: The name of the LC.',
               'COL time: The time point.',
               'COL time_unit: The unit of time used in time.',
               'COL percent_b: The fraction of solvent B in the LC output.',
               'COL flow: The solvent flow rate.',
               'COL flow_unit: The unit of flow used in flow']


class Inclusion(tables.IsDescription):
    """
    @brief HDF5 table definition for Inclusion list data table: data extracted from Xcalibur inclusion mass data
    """
    ms_mass = tables.Float32Col(pos=0)
    start = tables.Float32Col(pos=1)
    end = tables.Float32Col(pos=2)
    faims_cv = tables.Float32Col(pos=3)
    ms_col_energy = tables.Float32Col(pos=4)
    ms_charge = tables.Int32Col(pos=5)
    ms_intensity = tables.Float32Col(pos=6)
    ms2_mass = tables.Float32Col(pos=7)
    ms2_col_energy = tables.Float32Col(pos=8)
    name = tables.StringCol(100, pos=9)

Inclusion_attr = ['COL ms_mass: The m/z value to be included.',
                  'COL start: The start of the retention time window for inclusion.',
                  'COL end: The end of the retention time window for inclusion.',
                  'COL faims_cv: ?',
                  'COL ms_col_energy: The collision energy to be used.',
                  'COL ms_charge: The charge start of the ion.',
                  'COL ms_intensity: Intensity threshold for inclusion.',
                  'COL ms2_mass: The m/z value for a second fragmentation.',
                  'COL ms2_col_energy: The collision energy to be used for a second fragmentation.',
                  'COL name: Name of entitiy.']


class ExactiveInclusion(tables.IsDescription):
    """
    @brief HDF5 table definition for Exactive Inclusion list data table
    """
    ms_mass = tables.Float32Col(pos=0)
    start = tables.Float32Col(pos=1)
    end = tables.Float32Col(pos=2)
    polarity = tables.StringCol(3, pos=3)
    col_energy = tables.Float32Col(pos=4)
    charge = tables.Int32Col(pos=5)
    comment = tables.StringCol(100, pos=6)

ExactiveInclusion_attr = ['COL ms_mass: The m/z value to be included.',
                          'COL start: The start of the retention time window for inclusion.',
                          'COL end: The end of the retention time window for inclusion.',
                          'COL polarity: The polarity the instrument should be using to analyse the ion.',
                          'COL col_energy: The collision energy to be used.',
                          'COL charge: The charge start of the ion.',
                          'COL comment: Comments for this ion.']


class ExactiveLockMass(tables.IsDescription):
    """
    @brief HDF5 table definition for ExactiveLockMass data table
    """
    mass = tables.Float32Col(pos=0)
    start = tables.Float32Col(pos=1)
    end = tables.Float32Col(pos=2)
    polarity = tables.StringCol(3, pos=3)
    comment = tables.StringCol(100, pos=4)

ExactiveLockMass_attr = ['COL mass: The m/z value for the lock mass ion.',
                         'COL start: The start of the retention time window for inclusion.',
                         'COL end: The end of the retention time window for inclusion.',
                         'COL polarity: The polarity the instrument should be using to analyse the ion.',
                         'COL comment: Comments for this ion.']


class Reject(tables.IsDescription):
    """
    @brief HDF5 table definition for Reject masses table
    """
    mass = tables.Float32Col(pos=0)
    start = tables.Float32Col(pos=1)
    end = tables.Float32Col(pos=2)

Reject_attr = ['COL mass: The m/z value for the lock mass ion.',
               'COL start: The start of the retention time window for inclusion.',
               'COL end: The end of the retention time window for inclusion.']


class NeutralLoss(tables.IsDescription):
    """
    @brief closes the hdf5 file and creates indexes for farious data
    """
    mass = tables.Float32Col()

NeutralLoss_attr = ['COL mass: The m/z value for the neutral loss ion.']


class MSMSInclusion(tables.IsDescription):
    """
    @brief HDF5 table definition for MS/MS inclusion masses table
    """
    mode = tables.StringCol(1, pos=0)
    ms_mass = tables.Float32Col(pos=1)
    start = tables.Float32Col(pos=2)
    end = tables.Float32Col(pos=3)
    faims_cv = tables.Float32Col(pos=4)
    ms_col_energy = tables.Float32Col(pos=5)
    last_mass = tables.Float32Col(pos=6)
    name = tables.StringCol(100, pos=7)

MSMSInclusion_attr = ['COL mode: ?',
                      'COL ms_mass: The m/z value to be included.',
                      'COL start: The start of the retention time window for inclusion.',
                      'COL end: The end of the retention time window for inclusion.',
                      'COL faims_cv: ?',
                      'COL ms_col_energy: The collision energy to be used.',
                      'COL last_mass: ?',
                      'COL name: Name of entitiy.']


class Config(tables.IsDescription):
    """
    @brief HDF5 table definition for Configuration Parameters table
    """
    set = tables.StringCol(30, pos=0)
    parameter = tables.StringCol(60, pos=1)
    value = tables.StringCol(150, pos=2)

Config_attr = ['COL set: The config file group.',
               'COL parameter: The parameter name.',
               'COL value: The parameter value']


class Noise(tables.IsDescription):
    """
    @brief HDF5 table definition for noise data from MS and some MS/MS spectra
    """
    spec_id = tables.Int32Col(pos=0)
    mz = tables.Float32Col(pos=1)
    noise = tables.Float32Col(pos=2)
    baseline = tables.Float32Col(pos=3)

Noise_attr = ['COL spec_id: Spectrum identifier; corresponds to the Xcalibur spectrum number.',
              'COL mz: The m/z value for this measurement point.',
              'COL noise: The intensity cut off level at the m/z.  This is the best estimate of noise.',
              'COL baseline: ?']


class IsotopeLabel(tables.IsDescription):
    """
    @brief HDF5 table definition for reporter ion data
    """
    iso_id = tables.Int32Col(pos=0)
    name = tables.StringCol(10, pos=1)
    mz = tables.Float32Col(pos=2)
    intmz = tables.Int32Col(pos=3)
    error = tables.Float32Col(pos=4)
    method_id = tables.Int32Col(pos=5)
    amino = tables.StringCol(5, pos=6)

IsotopeLabel_attr = ['COL iso_id: The isotope labels id, should be unique amongst all methods.',
                     'COL name: The name for the individual isotope.',
                     'COL mz: The reporter ion m/z for the label for MS2 quantification, and the mass delta for MS1 quantification.',
                     'COL intmz: The mz rounded to the nearest integer.',
                     'COL error: The absolute m/z error to be applied to this label.',
                     'COL method_id: The unique identifer for the quantification method that this label is part of.',
                     'COL amino: The amino acid the label modifies.']


class Units(tables.IsDescription):
    """
    @brief HDF5 table definition for unit grouping of MS/MS scan events
    """
    unit = tables.Int32Col(pos=0)
    order = tables.Int32Col(pos=1)
    activation = tables.StringCol(8, pos=2)
    acttime = tables.Float32Col(pos=3)
    energy = tables.Float32Col(pos=4)
    isolation = tables.Float32Col(pos=5)
    lowmz = tables.Float32Col(pos=6)
    resolution = tables.Int32Col(pos=7)
    use = tables.StringCol(5, pos=8)
    scanevents = tables.StringCol(100, pos=9)
    colenergysteps = tables.Int32Col(pos=10)
    colenergywidth = tables.Float32Col(pos=11)

Units_attr = ['COL unit: The ID of the unit.',
              'COL order: The order number of the unit MS/MS event.',
              'COL activation: The fragmentation method used.',
              'COL acttime: The time allowed for collisions.',
              'COL energy: The energy used in the fragmentation.',
              'COL isolation: The isolation window setting.',
              'COL lowmz: Lowest m/z for the MS/MS event.',
              'COL resolution: Operating resolution of the instrument.',
              'COL use: How the spectrum data should be used I = identification only, Q = quantification only & IQ = quantification and identification.',
              'COL scanevents: The scan events from Xcalibur method linked to this .',
              'COL colenergysteps: The steps in collision energy if this option is being used.',
              'COL colenergywidth: Type range of collision energies used if multiple collision energies are being used.']


class Imports(tables.IsDescription):
    """
    @brief HDF5 table definition for mascot import tracking table
    """
    name = tables.StringCol(30, pos=1)
    date = tables.StringCol(30, pos=2)
    datfile = tables.StringCol(50, pos=3)
    status = tables.StringCol(10, pos=4)

Imports_attr = ['COL name: The name of the section in which the Mascot data is to be stored.',
                'COL date: Date and time of the Mascot import start.',
                'COL datfile: The name of the Mascot dat file to be imported.',
                'COL status: The status of the import.']


class Parameter(tables.IsDescription):
    """
    @brief HDF5 table definition for mascot parameters table
    """
    section = tables.StringCol(10, pos=0)
    parameter = tables.StringCol(30, pos=1)
    value = tables.StringCol(100, pos=2)

Parameter_attr = ['COL section: Section name from the MascotParser.cfg file.',
                  'COL parameter: The parameter name.',
                  'COL value: The parameter value']


class Mass(tables.IsDescription):
    """
    @brief HDF5 table definition for mascot masses table: Mascot elemental and amino acid masses
    """
    name = tables.StringCol(10, pos=0)
    mass = tables.Float32Col(pos=1)

Mass_attr = ['COL name: Name of the element or peptide terminal endings or single letter code of the amino acid.',
             'COL mass: Mass of the named element / peptide terminal endings / amino acid']


class Mod(tables.IsDescription):
    """
    @brief HDF5 table definition for mascot modifications table
    """
    id = tables.StringCol(1, pos=0)
    name = tables.StringCol(50, pos=1)
    modtype = tables.StringCol(10, pos=2)
    da_delta = tables.Float32Col(pos=3)
    amino = tables.StringCol(20, pos=4)
    neutralloss = tables.Float32Col(pos=5)
    nlmaster = tables.Float32Col(pos=6)
    relevant = tables.Int32Col(pos=7)
    composition = tables.StringCol(100, pos=8)

Mod_attr = ['COL id: Mascots identifier for the modification.  This is a number for variable modifications and an amino acid code for fixed modifications.',
            'COL name: The name of the modification.',
            'COL modtype: Whether the modification is fixed or variable.',
            'COL da_delta: The mass difference that the modification makes.',
            'COL amino: The amino acid specificity or terminal specificity.',
            'COL neutralloss: The masses of any neutral losses.',
            'COL nlmaster: The mass of the main neutral loss.',
            'COL relevant: Is the modification biologically relevant.',
            'COL composition: The elemental composition of the modification']


class Query(tables.IsDescription):
    """
    @brief HDF5 table definition for mascot query table: Mascot data from query and summary sections
    """
    query = tables.Int32Col(pos=0)
    spec_id = tables.Int32Col(pos=1)
    msms_id = tables.StringCol(10, pos=2)
    rt = tables.Float32Col(pos=3)
    prec_neutmass = tables.Float32Col(pos=4)
    prec_mz = tables.Float32Col(pos=5)
    prec_charge = tables.Int32Col(pos=6)
    matches = tables.Int32Col(pos=7)
    homology = tables.Float32Col(pos=8)
    numpeps = tables.Int32Col(pos=9)
    delta_seq = tables.Float32Col(pos=10)
    delta_mod = tables.Float32Col(pos=11)

Query_attr = ['COL query: The mascot assigned spectrum identifier.',
              'COL spec_id: Spectrum identifier; corresponds to the Xcalibur spectrum number.',
              'COL msms_id: spec_id reformatted in text format (F000000).',
              'COL rt: The retention time of the spectrum.',
              'COL prec_neutmass: The precursor m/z turned into zero charge state.',
              'COL prec_mz: The precursor m/z.',
              'COL prec_charge: The precursor charge state.',
              'COL matches: Number of peptides found by Mascot within the mass accuracy of the precursor.',
              'COL homology: The Mascot score required for a homology match.',
              'COL numpeps: A count of the number of peptides with valid sequences returned by Mascot for this spectrum.',
              'COL delta_seq: Difference in Mascot scores between the top Mascot sequence match and the next highest scoring sequence match. If delta_seq equals zero or if it is above the delta_seq_threshold its quant values may be used for protein quantification.',
              'COL delta_mod: Difference in Mascot scores between the best two modification variants of this sequence. See paper 21057138.']


class Peptide(tables.IsDescription):
    """
    @brief HDF5 table definition for mascot peptide table
    """
    query = tables.Int32Col(pos=0)
    pepno = tables.Int32Col(pos=1)
    sequence = tables.StringCol(255, pos=2)
    is_hook = tables.Int32Col(pos=3)
    useinprot = tables.Int32Col(pos=4)
    modsVariable = tables.StringCol(255, pos=5)
    modsFixed = tables.StringCol(255, pos=6)
    mass = tables.Float32Col(pos=7)
    da_delta = tables.Float32Col(pos=8)
    score = tables.Float32Col(pos=9)
    misscleave = tables.Int32Col(pos=10)
    numionsmatched = tables.Int32Col(pos=11)
    seriesfound = tables.StringCol(255, pos=12)
    peaks1 = tables.Int32Col(pos=13)
    peaks2 = tables.Int32Col(pos=14)
    peaks3 = tables.Int32Col(pos=15)

Peptide_attr = ['COL query: The Mascot assigned spectrum identifier.',
                'COL pepno: The Mascot assigned order of this peptide match to this query.',
                "COL sequence: Peptide's amino acid sequence.",
                'COL is_hook: Equals 1 if this is a hook peptide.',
                'COL useinprot: Equals 1 if the peptide is suitable to use in protein inference.',
                'COL modsVariable: Description of the Mascot identified variable modifications . e.g. 2 Oxidation (M); 1 TMT6plex (N-term).',
                'COL modsFixed: Description of the Mascot identified fixed modifications .  e.g. 1 Carbamidomethyl (C); 1 TMT6plex (K).',
                'COL mass: Calculated molecular mass of the peptide sequence.',
                'COL da_delta: Mass difference between the sequence mass (mw) and the neutral_mass in Da. Corresponds to Mascot delta.',
                'COL score: The Mascot score.',
                'COL misscleave: Number of sites on the peptide sequence where expected cleavage by the protease failed to occur.',
                'COL numionsmatched: The the number of ions from the sequence that match the spectrum data.',
                'COL seriesfound: Mascot seriesfound parameter.',
                'COL peaks1: Mascot peaks1 parameter.',
                'COL peaks2: Mascot peaks2 parameter.',
                'COL peaks3: Mascot peaks3 parameter.']


class ETpeptide(tables.IsDescription):
    """
    @brief HDF5 table definition for mascot error tolerant peptide table
    """
    query = tables.Int32Col(pos=0)
    pepno = tables.Int32Col(pos=1)
    sequence = tables.StringCol(255, pos=2)
    is_hook = tables.Int32Col(pos=3)
    useinprot = tables.Int32Col(pos=4)
    modsVariable = tables.StringCol(255, pos=5)
    modsFixed = tables.StringCol(255, pos=6)
    etModName = tables.StringCol(100, pos=7)
    etModDelta = tables.Float32Col(pos=8)
    mass = tables.Float32Col(pos=9)
    da_delta = tables.Float32Col(pos=10)
    score = tables.Float32Col(pos=11)
    misscleave = tables.Int32Col(pos=12)
    numionsmatched = tables.Int32Col(pos=13)
    seriesfound = tables.StringCol(255, pos=14)
    peaks1 = tables.Int32Col(pos=15)
    peaks2 = tables.Int32Col(pos=16)
    peaks3 = tables.Int32Col(pos=17)

ETpeptide_attr = ['COL query: SThe mascot assigned spectrum identifier.',
                  'COL pepno: The Mascot assigned order of this peptide match to this query.',
                  "COL sequence: Peptide's amino acid sequence.",
                  'COL is_hook: Equals 1 if this is a hook peptide.',
                  'COL useinprot: Equals 1 if the peptide is suitable to use in protein inference.',
                  'COL modsVariable: Description of the Mascot identified variable modifications . e.g. 2 Oxidation (M); 1 TMT6plex (N-term).',
                  'COL modsFixed: Description of the Mascot identified fixed modifications .  e.g. 1 Carbamidomethyl (C); 1 TMT6plex (K).',
                  'COL etModName: Name of error tollerant modification.',
                  'COL etModDelta: Mass difference of error tollerant modification).',
                  'COL mass: Calculated molecular mass of the peptide sequence.',
                  'COL da_delta: Mass difference between the sequence mass (mw) and the neutral_mass in Da. Corresponds to Mascot delta.',
                  'COL score: The Mascot score.',
                  'COL misscleave: Number of sites on the peptide sequence where expected cleavage by the protease failed to occur.',
                  'COL numionsmatched: The the number of ions from the sequence that match the spectrum data.',
                  'COL seriesfound: Mascot seriesfound parameter.',
                  'COL peaks1: Mascot peaks1 parameter.',
                  'COL peaks2: Mascot peaks2 parameter.',
                  'COL peaks3: Mascot peaks3 parameter.']


class Protein(tables.IsDescription):
    """
    @brief HDF5 table definition for mascot modifications table
    """
    accession = tables.StringCol(15, pos=0)
    mass = tables.Float32Col(pos=1)
    name = tables.StringCol(255, pos=2)

Protein_attr = ['COL accession: The protein identifier used by Mascot.',
                'COL mass: The mass of the protein sequence.',
                'COL name: The name of the protein.']


class Seq2Acc(tables.IsDescription):
    """
    @brief HDF5 table definition for links between peptide sequences and protein accessions
    """
    sequence = tables.StringCol(255, pos=0)
    accession = tables.StringCol(15, pos=1)
    hook = tables.Int32Col(pos=2)
    hookscore = tables.Float32Col(pos=3)
    pepscore = tables.Float32Col(pos=4)
    bestczrank = tables.Int32Col(pos=5)
    numpeps = tables.Int32Col(pos=6)
    start = tables.Int32Col(pos=7)
    end = tables.Int32Col(pos=8)

Seq2Acc_attr = ['COL sequence: The peptide amino acid sequence.',
                'COL accession: The protein accession.',
                'COL hook: Equals 1 if this peptide has been identified as a hook peptide.',
                'COL hookscore: The best score of any hook identifications of this peptide.',
                'COL pepscore: The best score of any identifications of this peptide.',
                'COL bestczrank: The best position (Mascot pepno) of any identifications of this peptide.',
                'COL numpeps: The number of times this peptide has been identified.',
                'COL start: The start location of the first occurrence of this peptide in the protein.',
                'COL end: The end location of the first occurrence of this peptide in the protein']


class Index(tables.IsDescription):
    """
    @brief HDF5 table definition for links between peptide sequences and protein accessions
    """
    section = tables.StringCol(20, pos=0)
    linenumber = tables.Int32Col(pos=1)

Index_attr = ['COL section: Mascot dat file section name.',
              'COL linenumber: The number of the line in the Mascot dat file where the section starts.']


class UMAminoAcids(tables.IsDescription):
    """
    @brief HDF5 table definition for unimod amino acid definitions
    """
    code = tables.StringCol(6, pos=0)
    name = tables.StringCol(20, pos=1)
    title = tables.StringCol(6, pos=2)
    monomass = tables.Float32Col(pos=3)
    avgmass = tables.Float32Col(pos=4)
    composition = tables.StringCol(100, pos=5)

UMAminoAcids_attr = ['COL code: Three letter amino acid code or terminus.',
                     'COL name: Full amino acid name or terminus.',
                     'COL title: One letter amino acid code or terminus (how the entity is used in Mascot.',
                     'COL monomass: The monoisotopic mass of the amino acid residue / terminus.',
                     'COL avgmass: The average mass of the amino acid residue / terminus.',
                     'COL composition: The elemental composition of the amino acid residue / terminus.']


class UMModifications(tables.IsDescription):
    """
    @brief HDF5 table definition for unimod modification definitions
    """
    modid = tables.Int32Col(pos=0)
    name = tables.StringCol(50, pos=1)
    title = tables.StringCol(50, pos=2)
    delta_mono_mass = tables.Float32Col(pos=3)
    delta_avge_mass = tables.Float32Col(pos=4)
    delta_comp = tables.StringCol(100, pos=5)

UMModifications_attr = ['COL modid: The Mascot modification identifier.',
                        'COL name: The name of the modification.',
                        'COL title: The short name as used by Mascot.',
                        'COL delta_mono_mass: The difference in monoisotopic mass that the modification makes.',
                        'COL delta_avge_mass: The difference in average mass that the modification makes.',
                        'COL delta_comp: The difference in elemental composition that the modification makes.']


class UMSpecificity(tables.IsDescription):
    """
    @brief HDF5 table definition for unimod modification specificity definitions
    """
    modid = tables.Int32Col(pos=0)
    group = tables.Int32Col(pos=1)
    site = tables.StringCol(20, pos=2)
    position = tables.StringCol(20, pos=3)
    classification = tables.StringCol(30, pos=4)
    hidden = tables.StringCol(5, pos=5)
    nl_mono_mass = tables.Float32Col(pos=6)
    nl_avge_mass = tables.Float32Col(pos=7)
    nl_comp = tables.StringCol(100, pos=8)

UMSpecificity_attr = ['COL modid: The Mascot modification identifier.',
                      'COL group: The Mascot group for the modification. This allows Mascot to consider multiple sites as one modification.',
                      'COL site: The amino acid or terminus where the modifcation can be found.',
                      'COL position: Location that the modification can occur within the peptide / protein.',
                      'COL classification: Biological classification of the potential source of this modification.',
                      'COL hidden: true / false if the modification is hidden and not to be used.',
                      'COL nl_mono_mass: The difference in monoisotopic mass that the neutral loss makes.',
                      'COL nl_avge_mass: The difference in average mass that the neutral loss makes.',
                      'COL nl_comp: The difference in elemental composition that the neutral loss makes.']


class UMElements(tables.IsDescription):
    """
    @brief HDF5 table definition for unimod element definitions
    """
    name = tables.StringCol(20, pos=0)
    title = tables.StringCol(5, pos=1)
    monomass = tables.Float32Col(pos=2)
    avgmass = tables.Float32Col(pos=3)

UMElements_attr = ['COL name: The name of the element.',
                   'COL title: The symbol used for the element.',
                   'COL monomass: The monoisotopic mass of the element.',
                   'COL avgmass: The average mass of the element.']


class Statistics(tables.IsDescription):
    """
    @brief HDF5 table definition for summary statistics
    """
    statistic = tables.StringCol(20, pos=0)
    value = tables.Float32Col(pos=1)

Statistics_attr = ['COL statistic: The name of the statistic.',
                   'COL value: The value of the statistic']


class ParameterWTypes(tables.IsDescription):
    """
    @brief HDF5 table definition for summary statistics
    """
    parameter = tables.StringCol(20, pos=0)
    type = tables.StringCol(5, pos=1)
    value = tables.StringCol(150, pos=2)

ParameterWTypes_attr = ['COL parameter: The name of the parameter.',
                        'COL type: The data type of the parameter (string, float etc.).',
                        'COL value: The value of the parameter.']


class Quan(tables.IsDescription):
    """
    @brief HDF5 table definition for quantification data table
    """
    spec_id = tables.Int32Col(pos=0)
    isolabel_id = tables.Int32Col(pos=1)
    area = tables.Float32Col(pos=2)
    corrected = tables.Float32Col(pos=3)
    inten = tables.Float32Col(pos=4)
    survey_id = tables.Int32Col(pos=5)
    mzdiff = tables.Float32Col(pos=6)
    ppm = tables.Float32Col(pos=7)
    coalescence = tables.Float32Col(pos=8)

Quan_attr = ['COL spec_id: Spectrum identifier; corresponds to the Xcalibur spectrum number.',
             'COL isolabel_id: The identifier of the isotope label.',
             'COL area: The area of the quantification ion.',
             'COL corrected: The intensity value after isotope interference correction.',
             'COL inten: The intensity value of the isotope label.',
             'COL survey_id: The MS1 spectrum prior to the MS/MS event.',
             'COL mzdiff: The m/z difference between the theoretical ion and the identified ion.',
             'COL ppm: mzdiff expressed in parts per million.',
             'COL coalescence: The measured degree of coalescence between adjacent label ions.  Scale between 0 no coalescnence and 1 fully coalesced.']


class CalMasses(tables.IsDescription):
    """
    @brief HDF5 table definition for calibration masses table
    """
    query = tables.Int32Col(pos=0)
    pepno = tables.Int32Col(pos=1)
    spec_id = tables.Int32Col(pos=2)
    rt = tables.Float32Col(pos=3)
    sequence = tables.StringCol(255, pos=4)
    modsVariable = tables.StringCol(255, pos=5)
    modsFixed = tables.StringCol(255, pos=6)
    mass = tables.Float32Col(pos=7)
    da_delta = tables.Float32Col(pos=8)
    score = tables.Float32Col(pos=9)
    prec_neutmass = tables.Float32Col(pos=10)
    prec_mz = tables.Float32Col(pos=11)
    prec_charge = tables.Int32Col(pos=12)
    hitType = tables.StringCol(255, pos=13)

CalMasses_attr = ['COL query: SThe mascot assigned spectrum identifier.',
                  'COL pepno: The Mascot assigned order of this peptide match to this query.',
                  'COL spec_id: Spectrum identifier; corresponds to the Xcalibur spectrum number.',
                  'COL rt: The retention time of the spectrum.',
                  "COL sequence: Peptide's amino acid sequence.",
                  'COL modsVariable: Description of the Mascot identified variable modifications . e.g. 2 Oxidation (M); 1 TMT6plex (N-term).',
                  'COL modsFixed: Description of the Mascot identified fixed modifications .  e.g. 1 Carbamidomethyl (C); 1 TMT6plex (K).',
                  'COL mass: Calculated molecular mass of the peptide sequence.',
                  'COL da_delta: Mass difference between the sequence mass (mw) and the neutral_mass in Da. Corresponds to Mascot delta.',
                  'COL score: The Mascot score.',
                  'COL prec_neutmass: The precursor m/z turned into zero charge state.',
                  'COL prec_mz: The precursor m/z.',
                  'COL prec_charge: The precursor charge state.',
                  'COL hitType: Type of match.']


class ResultSample(tables.IsDescription):
    '''
    @brief HDF5 table definition for teh stand alone results Sample table
    '''

    sample_id = tables.Int32Col(pos=0)
    source_file = tables.StringCol(150, pos=1)
    acquisition_time = tables.Float32Col(pos=2)
    acquired_spectra = tables.Int32Col(pos=3)
    mascot_matched_spectra = tables.Int32Col(pos=4)
    spectra_in_qc_proteins = tables.Int32Col(pos=5)
    quantified_spectra = tables.Int32Col(pos=6)
    mean_precursor_ion_accuracy = tables.Float32Col(pos=7)
    sd_precursor_ion_accuracy = tables.Float32Col(pos=8)
    mean_reporter_ion_accuracy = tables.Float32Col(pos=9)
    sd_reporter_ion_accuracy = tables.Float32Col(pos=10)


ResultSample_attr = ['COL sample_id: Generated identifier of the sample where the spectrum was measured. Master entry.',
                     'COL source_file: The file name of the original Xcalibur .raw file name from which the data were acquired.',
                     'COL acquisition_time: The total acquisition time of the Xcalibur analysis ',
                     'COL acquired_spectra: The number of MS/MS events acquired in the Xcalibur raw file.',
                     'COL mascot_matched_spectra: The number of spectra  for which Mascot found a peptide match.',
                     'COL spectra_in_qc_proteins: The number of spectra with peptides linked to proteins fulfilling all QC criteria.',
                     'COL quantified_spectra: The number of spectra fulfilling all quantification-specific QC criteria linked to validated proteins.',
                     'COL mean_precursor_ion_accuracy: The precursor ion accuracy is calculated as: (measured precursor ion m/z - theoretical peptide m/z) / theoretical peptide m/z, expressed in ppm.  The mean is calculated from all spectra linked to accepted proteins.',
                     'COL sd_precursor_ion_accuracy: The standard deviation of the data calculated for the mean_precursor_ion_accuracy, expressed in Th.',
                     'COL mean_reporter_ion_accuracy: The reporter ion accuracy is calculated as: measured m/z - theoretical m/z, expressed in Th.  The mean is calculated from all detected reporter ions from all spectra with a linked to a valid protein.',
                     'COL sd_reporter_ion_accuracy: The standard deviation of the data calculated for the mean_reporter_ion_accuracy']


class ResultSpectrum(tables.IsDescription):
    '''
    @brief HDF5 table definition for teh stand alone results Spectrum table
    '''

    spectrum_id = tables.Int32Col(pos=0)
    sample_id = tables.Int32Col(pos=1)
    msms_id = tables.Int32Col(pos=2)
    query = tables.Int32Col(pos=3)
    neutral_mass = tables.Float32Col(pos=4)
    parent_ion = tables.Float32Col(pos=5)
    precursor_mz = tables.Float32Col(pos=6)
    s2i = tables.Float32Col(pos=7)
    p2t = tables.Float32Col(pos=8)
    survey_id = tables.Int32Col(pos=9)
    start_time = tables.Float32Col(pos=10)
    peak_rt = tables.Float32Col(pos=11)
    peak_intensity = tables.Float32Col(pos=12)
    peak_fwhm = tables.Float32Col(pos=13)
    sum_reporter_ions = tables.Float32Col(pos=14)
    quant_cancelled = tables.Int32Col(pos=15)
    charge_state = tables.Int32Col(pos=17)


ResultSpectrum_attr = ['COL spectrum_id: Assigned identifier for this spectrum.  Is unique across all merged datasets in a single analysis.  Master entry.',
                       'COL sample_id: Generated identifier of the sample where the spectrum was measured.',
                       'COL msms_id: Identifier of MS/MS spectrum from the raw data file. Non unique across merged datasets.',
                       'COL query: The spectrum identifier as assigned by Mascot.',
                       'COL neutral_mass: The neutral mass calculated from the precursor_mz by Mascot. Corresponds to neutral_mass in Mascot.',
                       'COL parent_ion: The actual m/z setting used by the instrument for the isolation of the precursor ion.',
                       'COL precursor_mz: Observed m/z of precursor ion. Corresponds to "observed" in Mascot .dat file.',
                       'COL s2i: Signal-to-interference value for the precursor ion of this spectrum',
                       'COL p2t: Intensity- to -noise value for the precursor ion of this spectrum.  (p2t = precursor-to-threshold)',
                       'COL survey_id: Identifier of the survey spectrum from the .raw file. Non unique across merged datasets.',
                       'COL start_time: Retention time  for the MS/MS event',
                       'COL peak_rt: Retention time of the apex of the chromatographic peak identified linked to this spectrum.',
                       'COL peak_intensity: Maximum intensity of the XIC peak identified from this spectrum.',
                       'COL peak_fwhm: The full width at half maximum intensity of the chromatographic peak identified from this spectrum.',
                       'COL sum_reporter_ions: Sum of all reporter ion intensities in this spectrum.',
                       'COL quant_cancelled: Is the quantification of this spectrum cancelled',
                       'COL charge_state: Charge state of the precursor ion.']


class ResultSpecQuant(tables.IsDescription):
    '''
    @brief HDF5 table definition for the stand alone results SpecQuant table
    '''

    spectrum_id = tables.Int32Col(pos=0)
    isotopelabel_id = tables.Int32Col(pos=1)
    protein_group_no = tables.Int32Col(pos=2)
    quant_raw = tables.Float32Col(pos=3)
    quant_isocorrected = tables.Float32Col(pos=4)
    quant_allcorrected = tables.Float32Col(pos=5)
    score = tables.Float32Col(pos=6)
    delta_seq = tables.Float32Col(pos=7)
    s2i = tables.Float32Col(pos=8)
    p2t = tables.Float32Col(pos=9)
    in_quantification_of_protein = tables.Int32Col(pos=10)
    peptide_length = tables.Int32Col(pos=11)
    is_unique = tables.Int32Col(pos=12)
    fdr_at_score = tables.Float32Col(pos=13)


ResultSpecQuant_attr = ['COL spectrum_id: Assigned identifier for this spectrum.  Is unique across all merged datasets in a single analysis. Corresponds to master entry in spectrum table.',
                        'COL isotopelabel_id: Isotope label identifier.  These are defined in the QuantMethod.cfg file and should be unique across all quantification methods.',
                        'COL protein_group_no: Generated protein number to which the peptide sequence is matched. Corresponds to protein entry in proteins output.',
                        'COL quant_raw: Unprocessed quantification value.  In MS2 quantification this is the intensity of the reporter ion.',
                        'COL quant_isocorrected: Quantification value following subtraction of isotope interference from adjacent label(s).',
                        'COL quant_allcorrected: Quantification value after all corrections have been applied.',
                        'COL score: Mascot score',
                        'COL delta_seq: Difference in Mascot scores between this Mascot sequence match and the next highest scoring sequence match. If delta_seq equals zero or if it is above the delta_seq_threshold its quant values may be used for protein quantification.',
                        "COL s2i: Signal-to-interference value for the precursor ion of this spectrum. If the s2i value is above s2ithreshold the peptide's quantification values may be used for protein quantification.",
                        'COL p2t: Intensity-to-noise value for the precursor ion of this spectrum. If the p2t value is above p2tthreshold the peptides quantification values may be used for protein quantification. (p2t = precursor2threshold).',
                        "COL in_quantification_of_protein: Equal to 1 if peptide's reporter ions are used for protein quantification.",
                        'COL fdr_at_score: Calculated peptide false discovery rate for this Mascot score. This may be used as a threshold to exclude peptide from use in protein quantification and from protein inference.',
                        'COL peptide_length: Length of the peptide sequence.',
                        'COL is_unique: Peptide sequence is found uniquely in given protein group. This is required for its quant values to be used in protein quantification.']


class ResultPeptide(tables.IsDescription):
    '''
    @brief HDF5 table definition for teh stand alone results Peptide table
    '''

    peptide_id = tables.Int32Col(pos=0)
    spectrum_id = tables.Int32Col(pos=1)
    protein_group_no = tables.Int32Col(pos=2)
    peptide = tables.StringCol(100, pos=3)
    variable_modstring = tables.StringCol(100, pos=4)
    fixed_modstring = tables.StringCol(100, pos=5)
    positional_modstring = tables.StringCol(350, pos=6)
    score = tables.Float32Col(pos=7)
    rank = tables.Int32Col(pos=8)
    mw = tables.Float32Col(pos=9)
    da_delta = tables.Float32Col(pos=10)
    ppm_error = tables.Float32Col(pos=11)
    missed_cleavage_sites = tables.Int32Col(pos=12)
    delta_seq = tables.Float32Col(pos=13)
    delta_mod = tables.Float32Col(pos=14)
    is_hook = tables.Int32Col(pos=15)
    is_duplicate = tables.Int32Col(pos=16)
    is_first_use_of_sequence = tables.Int32Col(pos=17)
    is_unique = tables.Int32Col(pos=18)
    is_reverse_hit = tables.Int32Col(pos=19)
    is_quantified = tables.Int32Col(pos=20)
    failed_fdr_filter = tables.Int32Col(pos=21)
    fdr_at_score = tables.Float32Col(pos=22)
    in_protein_inference = tables.Int32Col(pos=23)
    seq_start = tables.Int32Col(pos=24)
    seq_end = tables.Int32Col(pos=25)


ResultPeptide_attr = ['COL peptide_id: Generated unique identifier for this peptide record.',
                      'COL spectrum_id: Assigned identifier for this spectrum.  Is unique across all merged datasets in a single analysis. Corresponds to master entry in spectrum table.',
                      'COL protein_group_no: Generated protein number to which the peptide sequence is matched. Corresponds to master entry in proteinhit table.',
                      "COL peptide: Peptide's amino acid sequence. The length of the peptide sequence is used in peptide QC and to exclude it from use in protein quantification. This is given as sequence in the output.",
                      'COL variable_modstring: Description of the Mascot identified variable modifications . e.g. 2 Oxidation (M); 1 TMT6plex (N-term).',
                      'COL fixed_modstring: Description of the Mascot identified fixed modifications .  e.g. 1 Carbamidomethyl (C); 1 TMT6plex (K).',
                      'COL positional_modstring: Description of all modifications with their position in the sequence .  e.g. TMT6plex:0; Carbamidomethyl:2; Oxidation:3; Oxidation:10; TMT6plex:12.',
                      'COL score: Score assigned to spectrum-to-protein match by Mascot. This may be used as a threshold to exclude peptide from being used for protein quantification.',
                      'COL rank: Rank of spectrum to protein match within Mascots list of ten suggestions. Only rank1 peptides are quantified.',
                      'COL mw: Calculated molecular mass of the peptide sequence corresponds to mrcalc in Mascot .dat file.',
                      'COL da_delta: Mass difference between the sequence mass (mw) and the neutral_mass in Da. Corresponds to Mascot delta.',
                      'COL ppm_error: da_error expressed in ppm.',
                      'COL missed_cleavage_sites: Number of sites on the peptide sequence where expected cleavage by the protease failed to occur.',
                      'COL delta_seq: Difference in Mascot scores between this Mascot sequence match and the next highest scoring sequence match.',
                      'COL delta_mod: Difference in Mascot scores between the best two modification variants of this sequence. See ref - pubmed #21057138.',
                      'COL is_hook: Equals 1 if this is a hook peptide.',
                      'COL is_duplicate: When multiple identifications are made to the same sequence the highest scoring identification has is_duplicate =0. All other identifications have is_duplicate = 1.',
                      'COL is_first_use_of_sequence: Equals 1  for the first occurrence of this peptide sequence across all protein groups when they are ordered by descending total_score.',
                      'COL is_unique: Peptide sequence is found uniquely in given protein group. This is required for its quant values to be used in protein quantification.',
                      'COL is_reverse_hit: Equals 1 if peptide sequence maps only to a protein from the reverse database.',
                      'COL is_quantified: Equals 1 if peptide links to at least one quantification event (reporter ion)',
                      "COL failed_fdr_filter: Equals 1 if peptide's FDR is above the threshold given at runtime. This means it is not used for quantification or for protein inference during the aggregation / re-assignment step.",
                      'COL fdr_at_score: Calculated peptide false discovery rate for this score. This may be used as a threshold to exclude peptide from use in protein quantification and from protein inference.',
                      'COL in_protein_inference: Equals 1 if peptide sequence is used in assignment of protein groups.',
                      'COL seq_start: Peptide start position on protein sequence given by first identifying accession in protein group.',
                      'COL seq_end: Peptide end position on protein sequence given by first identifying accession in protein group.']


class ResultProteinHit(tables.IsDescription):
    '''
    @brief HDF5 table definition for teh stand alone results ProteinHit table
    '''

    protein_group_no = tables.Int32Col(pos=0)
    protein_id = tables.StringCol(100, pos=1)
    description = tables.StringCol(500, pos=2)
    gene_name = tables.StringCol(50, pos=3)
    mw = tables.Float32Col(pos=4)
    ssm = tables.Int32Col(pos=5)
    total_score = tables.Float32Col(pos=6)
    hssm = tables.Int32Col(pos=7)
    hookscore = tables.Float32Col(pos=8)
    upm = tables.Int32Col(pos=9)
    is_reverse_hit = tables.Int32Col(pos=10)
    protein_fdr = tables.Float32Col(pos=11)
    max_score = tables.Float32Col(pos=12)


ResultProteinHit_attr = ['COL protein_group_no: Generated protein number for this protein group. Master entry.',
                         'COL protein_id: Identifying accession(s) from protein database. In cases where all identified peptides match to more than one protein sequence, all possible accessions are given here. The order in which the accessions are presented determines the order of the corresponding information in the description and mw columns. This protein_id can be used to match proteins across files from multiple experiments (as opposed to protein_group_no).',
                         'COL description: Description line from .fasta file used in Mascot search. More than one entry possible.',
                         'COL gene_name: For data searched with Uniprot fasta files only. This corresponds to the GN value given in the description line of the fasta file. For multiple peptide match entries this value is the minimal set of gene names associated with the entry. This may be a more appropriate identifier used for matching proteins between files from multiple experiments.',
                         'COL mw: Molecular weight of full protein sequences. More than one entry possible see definition of protein_id.',
                         'COL ssm: Number of spectrum-to-sequence matches [peptides] in protein group (ssm=spectrum-sequence matches).',
                         'COL total_score: Sum of mascot scores for all peptides in protein group greater than peptide FDR threshold. Note that only highest Mascot score is counted for all versions of peptide sequence in that group',
                         'COL hssm: Number of spectra matched to hook peptides in protein group. See definition of hook peptide above. hssm = hook sequence-spectra matches.',
                         'COL hookscore: Sum of the Mascot scores for the matched Hook peptides',
                         'COL upm: Number of peptide sequences in protein group that are unique to it (upm=unique peptide matches).',
                         'COL is_reverse_hit: is the match a false positive match',
                         'COL protein_fdr: False discovery rate for protein based on max_score.',
                         'COL max_score: Maximum Mascot score of all peptides in protein group.']


class ResultStatistics(tables.IsDescription):
    '''
    @brief HDF5 table definition for the stand alone results statistics table
    '''

    statistic = tables.StringCol(100, pos=0)
    value = tables.StringCol(200, pos=1)


ResultStatistics_attr = ['COL statistic: Name of summary statistic',
                         'COL value: Value of summary statistic']


class ResultProteinQuant(tables.IsDescription):
    '''
    @brief HDF5 table definition for teh stand alone results ProteinQuant table
    '''

    protein_group_no = tables.Int32Col(pos=0)
    isotopelabel_id = tables.Int32Col(pos=1)
    reference_label = tables.StringCol(50, pos=2)
    protein_fold_change = tables.Float32Col(pos=3)
    lower_confidence_level = tables.Float32Col(pos=4)
    upper_confidence_level = tables.Float32Col(pos=5)
    sum_quant_signal = tables.Float32Col(pos=6)
    qssm = tables.Int32Col(pos=7)
    qupm = tables.Int32Col(pos=8)


ResultProteinQuant_attr = ['COL protein_group_no: Generated protein number to which the peptide sequence is matched. Corresponds to master entry in proteinhit table.',
                           'COL isotopelabel_id: Isotope label identifier.  These are defined in the QuantMethod.cfg file and should be unique across all quantification methods.',
                           'COL reference_label: Name of reporter ion label corresponding to the control on which relative protein quantification is based.',
                           'COL protein_fold_change: Protein fold change compared to reference.',
                           'COL lower_confidence_level: Lower bound of confidence level interval for the protein fold change.',
                           'COL upper_confidence_level: Upper bound of confidence level interval for the protein fold change.',
                           'COL sum_quant_signal: Sum of all reporter ion signals of valid quantifiable peptides in protein group following all relevant corrections.',
                           'COL qssm: Number of spectrum to sequence matches [peptides] in protein group with reporter ions used in protein quantification (qssm =quantfied spectrum-sequence matches).',
                           'COL qupm: Number of unique peptide sequences in protein group with reporter ions used in protein quantification (qupm=quantified unique peptide matches).']


class FDRData(tables.IsDescription):
    data_type = tables.StringCol(50, pos=0)
    score = tables.Float32Col(pos=1)
    reverse_hits = tables.Int32Col(pos=2)
    forward_hits = tables.Int32Col(pos=3)
    global_fdr = tables.Float32Col(pos=4)
    local_fdr = tables.Float32Col(pos=5)
    true_spectra = tables.Int32Col(pos=6)


FDRData_attr = ['COL data_type: Describes whether the data is a peptide or protein based FDR',
                'COL score: Mascot score',
                'COL reverse_hits: Number of reverse (false) hits at Mascot score',
                'COL forward_hits: Number of forward (true) hits at Mascot score',
                'COL global_fdr: FDR for Mascot score >= score',
                'COL local_fdr: FDR for Mascot score = score',
                'COL true_spectra: Number of spectra with Mascot score >= score']


class ResultConfigParameters(tables.IsDescription):
    """
    @brief HDF5 table definition for postMascot processing Parameters table
    """
    analysis = tables.StringCol(30, pos=0)
    application = tables.StringCol(30, pos=1)
    section = tables.StringCol(30, pos=2)
    parameter = tables.StringCol(60, pos=3)
    value = tables.StringCol(100, pos=4)


ResultConfigParameters_attr = ['COL analysis: Name of the part of the analysis that is calling the application.',
                               'COL application: Name of the application using the parameter',
                               'COL section: Section of the config file containing the parameter',
                               'COL parameter: Parameter name', 'COL value: Parameter value']
