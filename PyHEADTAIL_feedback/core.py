import numpy as np
import types
import copy
version = '0.2.1.1'

"""
@author Jani Komppula
@date: 11/10/2017

This file contains the core functions and variables for the signal processing
framework.

The basic concept is that a signal compatible with the framework
is generated by using interfaces (e.g. see feedback.py for PyHEADTAIL) which
utilize the basic functions and tools from the core (this file). A signal
from the interfaces is processed by passing it through a list of signal
processors, which represents, for example, a model of a transverse feedback
system. The signal processors in the list represent elementary analog and
digital signal processing steps, e.g. from the pickup plate to the kicker.
The signal processing model of the system itself is independent of the
interfaces, i.e. it can be used with any signal source (e.g. PyHEADTAIL or
signal tools for testing).

Both interfaces and signal processors can be developed separately and
dynamically without interfering with each other. This file contains an
example of a signal processor. More signal processors can be found from
the folder 'processors'

This file has been divided into three sections. The sections "SIGNALS AND
PARAMETERS" and "SIGNAL PROCESSORS" contain specifications, generators and
examples for signals and signal processors. The core functions for processing
signals with the signal processors, which can be used in programming
interfaces and signal processors are presented in the section "TOOLS".
"""


"""
## SIGNALS AND PARAMETERS
=========================
The definition of signal depends on context. For example, in one context, a
signal might be a continuous, time-varying, quantity, i.e. an analog signal. In
another context, it might be a list of numbers which represent values of
certain moments of equally spaced time, i.e. a digital signal. There are also
situations, where signals are more complex.

Unfortunately, the framework should work with all signals in all the mentioned
contexts. Because it is very unpractical if not almost impossible to develop
all signal processors to work with the all types of signals, different
definitions for a signal can be used in the framework.

The basic definition of a signal is that it consists of discrete numbers
in time domain. More specifically, a signal is a numpy array of
floating point numbers in time domain, i.e.
"""


def Signal(signal=[]):
    """ Returns a prototype for a signal."""
    return np.array(signal)


"""
Each number in the array corresponds to a signal value in the
specific time interval, *bin*. By default, the value of the bin
is a time average of the signal over the bin, but this is not
guaranteed because normalization of the signal processors depends
on the studied case.

In order to simplify development of signal processors, signals are
categorized into three classes basing on the assumptions which can be made
from a signal. Thus, signal processors can be specified to receive and
transmit signals in the specific classes. Due to the hierarchy of the signal
classes, a signal processor designed for lower class signals is also able to
process signals from higher classes.

    ### Class 0 signals
    -------------------
    There are no limitations for Class 0 signals, i.e. bin spacing and bin
    length might vary randomly. If the signal can be divided into segments,
    each segment must have an equal number of bins with equal bin spacing and
    bin lengths.

    Class 0 signal gives a large freedom to use any kind of signal as an input
    for the signal processors. Particularly it means that a single array of
    the slice values from multiple bunches can be used directly as a signal.

    ### Class 1 signals
    -------------------
    In this class, it is assumed that signal might be divided into equal
    length sequences which are separated by empty spaces. Bin spacing and
    width must be constant and equal in each segment.

    In practice, this means that signals from each bunch have an equal number
    of equally spaced slices/samples.

    ### Class 2 signals
    -------------------
    The signal is equally spaced and continuous in time.

    In practice, this means that the signal is continuously sliced/sampled over
    all bunches including empty spaces between bunches. This also limits the
    slicing/sampling rate to be a fraction of the bunch spacing in the case of
    multi bunch simulations.

The signal itself does not contain any information about what is the signal
class or how the bins are located in the physical space. Thus, this
information is given in parallel to the signal to the signal processors
by using a dictionary *parameters*.

The standard (minimal) prototype for the parameter dictionary is following,
but it is also allowed to add additional information to be carried between
the signal processors:
"""

def Parameters(signal_class=0, bin_edges=np.array([]), n_segments=0,
               n_bins_per_segment=0, segment_ref_points=np.array([]),
               previous_parameters = [], location=0, beta=1.):
    """
    Returns a prototype for signal parameters.

    Parameters
    ----------
    class : int
        A signal class
    bin_edges : NumPy array
        A 2D numpy array, which is equal length to the signal. Each row
        includes two floating point numbers, the edge positions of
        the bin in the physical space (time [s]).
    n_segments : int
        A number of equal length and equally binned segments where to
        the signal can be divided
    n_bins_per_segment : int
        A number of bins per segment. `len(bin_edges)/n_segments`
    segment_ref_points : NumPy array
        A numpy array of the reference point for the segments
    previous_parameters : array
        A list of Parameter objects, which tracks how the samping is changed
        during the signal processing
    location : float
        A location of the signal in betatron phase.
    beta : float
        A vale of beta function in the source of the signal. Value 1
        is neutral for signal processing
    """
    return {'class': signal_class,
            'bin_edges': bin_edges,
            'n_segments': n_segments,
            'n_bins_per_segment': n_bins_per_segment,
            'segment_ref_points': segment_ref_points,
            'previous_parameters': previous_parameters,
            'location': location,
            'beta': beta
            }


"""
In prinicple, the framework can also be extended to work
in other domains (e.g. f, s or z, symbolic calculations or even influences
from functional programmin), if the definition of signal is extended and
specific processors for the signal conversion between the different domains are
programmed.

"""

"""
## SIGNAL PROCESSORS
====================
A signal processor is a Python object which processes/modifies signals. The
signal is processed in the method process(parameters, signal, *args, **kwargs),
which takes arguments *parameters* and *signal* and returns (possibly) modified
versions of them.

The constructor of a signal processor must include the following lines:
    * Signal classes for the input and output signals
        self.signal_classes = (int, int)
    * A list of possible extensions used in the processors (an empty list by
      default):
        self.extensions = []    
    * Default macros, which help debugging and assist the future development
        self._macros = [] + default_macros(self, 'ProcessorName', **kwargs)

An example code for a minimal signal processor, which only bypasses a signal is
following:
"""
class Bypass(object):
    def __init__(self, **kwargs):
        self.signal_classes = (0, 0)
        self.extensions = []
        # a list of macros
        self._macros = [] + default_macros(self, 'Bypass', **kwargs)

    def process(self, parameters, signal, *args, **kwargs):

        return parameters, signal

"""
The framework supports extensions to the minimial layout. For example,
extensions can be used to provide extra simulation data to signal processors
(bunch extension), implement more complex signal transfer between simulation
objects (registers and combiners) or provide extra data for debugging and
data visualization. Names of the extensions supported by the core are listed
on the variable extensions.

The following extensions are supported by the core of the
framework.

    ### Bunch extension
    -------------------
    The bunch extensions allows use of additional slice set data from PyHEADTAIL
    in the processor. Because calculations of the statistical variables of the
    slice sets require a lot of computing power, the names of required
    variables are listed on the variable required_variables.

    The slice set data can be found from kwargs['slice_sets'], which is
    a list of slice set objects (emulations) of the simulated bunches. Note
    that it is not quaranteed that the bin set of the signal corresponds to
    the bin set used in the slice set. This can be checked by checking a
    number of items in the element 'previous_parameters' of the input
    parameters.

"""
class MinimalChargeWeighter(object):
    def __init__(self, **kwargs):
        # Signal classes for incoming and outgoing signals
        self.signal_classes = (0, 0)

        # A list of extensions supported by the processor
        self.extensions = ['bunch']

        # A list of PYHEADTAIL slice set variables required by the processor
        self.required_variables = ['n_macroparticles_per_slice']

    def process(self, parameters, signal, *args, **kwargs):

        slice_sets = kwargs['slice_sets']

        output_signal = np.copy(signal)

        for i, slice_set in enumerate(slice_sets):
            n_macroparticles = np.sum(slice_set.n_macroparticles_per_slice)

            j_from = i * parameters['n_bins_per_segment']
            j_to = (i + 1) * parameters['n_bins_per_segment']
            output_signal[j_from:j_to] *= slice_set.n_macroparticles_per_slice/n_macroparticles

        # The signal or the prameters could be modified here
        return parameters, output_signal


"""
    ### Register and combiner extensions
    ----------------------
    Register and combiner processors are designed to save and combine data
    from multiple turns. The basic principle of a register is that it is an
    iterable object, i.e. data from the register can be riden by iterating it
    in a for loop. As a part of the signal processor list, register by passes
    the signal without modifications. A combiner uses a register as a signal
    source and returns a combined signal, i.e. it reads the values from
    register(s) and combines signals by applying possible betatron phase
    advance correction algorithms.

    Details of these processors can be found from the
    file processors/register.py.
"""


"""
## TOOLS
========
"""

def process(parameters, signal, processors, **kwargs):
    """
    Processes a signal through the given signal processors

    Parameters
    ----------
    parameters : dict
        A standardized dict of the additional parameters describing the signal
    signal : NumPy array
        The signal
    processors : list
        A list of signal processors.
    **kwargs : -
        Other arguments which will be passed to the signal processors

    Returns
    -------
    dict
        Possibly modified dict of the signal parameters
    NumPy array
        The processed signal
    """

    for processor in processors:
        parameters, signal = processor.process(parameters, signal, **kwargs)
#        if signal is None:
#            print 'None signal!'
#            break

    return parameters, signal


def bin_widths(bin_edges):
    return (bin_edges[:, 1]-bin_edges[:, 0])


def bin_mids(bin_edges):
    return (bin_edges[:, 0]+bin_edges[:, 1])/2.


def bin_edges_to_z_bins(bin_edges):
    return np.append(bin_edges[:, 0], bin_edges[-1, 1])


def z_bins_to_bin_edges(z_bins):
    return np.transpose(np.array([z_bins[:-1], z_bins[1:]]))


def append_bin_edges(bin_edges_1, bin_edges_2):
    return np.concatenate((bin_edges_1, bin_edges_2), axis=0)


def get_processor_extensions(processors, external_extensions=None):
    """
    A function, which checks available extensions from the processors.

    Parameters
    ----------
    processors : list
        A list of signal processors.
    external_extensions : list
        A list of external extensions, which will be added to the list

    Returns
    -------
    list
        A list of found extensions
    """

    if external_extensions is None:
        available_extensions = []
    else:
        available_extensions = external_extensions

    for processor in processors:
        if processor.extensions is not None:
            available_extensions.extend(processor.extensions)

    available_extensions = list(set(available_extensions))

    return available_extensions

"""
### Extension specific functions
================================
"""
def get_processor_variables(processors, required_variables=None):
    """
    A function which checks the required PyHEADTAIL slice set variables
    from the signal processors.

    Parameters
    ----------
    processors : list
        A list of signal processors.
    external_variables : list
        A list of external variables, which will be added to the list

    Returns
    -------
    list
        A list of found statistical variables
    """

    if required_variables is None:
        required_variables = []

    for processor in processors:
        if 'bunch' in processor.extensions:
            required_variables.extend(processor.required_variables)

    required_variables = list(set(required_variables))

    if 'z_bins' in required_variables:
        required_variables.remove('z_bins')

    return required_variables


"""
### MACROS
========
"""
def default_macros(obj, label=None, **kwargs):
    func_list = []

    func_list = func_list + debug_macro(obj, label=label, **kwargs)
    func_list = func_list + label_macro(obj, label=label, **kwargs)
    func_list = func_list + init_vatiables_macro(obj, **kwargs)

    return func_list


def label_macro(obj, label=None, **kwargs):
    setattr(obj, 'label', label)
    return []


def init_vatiables_macro(obj, **kwargs):
    setattr(obj, 'time_scale', 0)
    return []


def debug_macro(obj, **kwargs):
    """
    A debug macro.
    
    If input parameter debug = True is given to the signal
    processor, the input and output parameters and signals are stored
    to the signal processor.

    Parameters
    ----------
    target_object : object
        A object which is operated (virtually always self)
    label : string
        A name of the signal processor

    Returns
    -------
    object
        A macro function, which is run in the process(...) method in the debug
        mode is turn on.
    """
    def decorated_process(self, parameters, signal, *args, **kwargs):
        input_parameters = parameters
        input_signal = np.copy(signal)

        output_parameters, output_signal = self.process_org(parameters, signal,
                                                       *args, **kwargs)

        for macro in self._macros:
            macro(self, input_parameters, input_signal, output_parameters,
                  output_signal, *args, **kwargs)

        return output_parameters, output_signal

    def store_data(target_object, input_parameters, input_signal,
                   output_parameters, output_signal, *args, **kwargs):
        if target_object.debug:
            target_object.input_parameters = copy.copy(input_parameters)
            target_object.input_signal = np.copy(input_signal)
            target_object.output_parameters = copy.copy(output_parameters)
            target_object.output_signal = np.copy(output_signal)

    if 'debug' in kwargs:
        obj.extensions.append('debug')

        setattr(obj, 'debug', kwargs['debug'])
        setattr(obj, 'input_parameters', None)
        setattr(obj, 'input_signal', None)
        setattr(obj, 'output_parameters', None)
        setattr(obj, 'output_signal', None)

        obj.process_org = obj.process
        obj.process = types.MethodType(decorated_process, obj)

        return [store_data]
    return []
